# Load cases on the pure path (CaseLoads, case_loads, case_Pf).
# C0 exit criteria: the table reproduces the stateful per-case load vectors
# exactly, the design-dependent fixed-end forces of every case re-evaluate to
# the cached ones at the model's own state and to the single-case
# `assemble_Pf` at any other state, and the cache is left untouched.

using LinearAlgebra, SparseArrays, Zygote, ForwardDiff, FiniteDifferences, ChainRulesCore

# a slab of two quads and four triangles on an edge beam, a haunched
# (VariableElement) beam and a column; three cases mixing every load kind.
# `variable = false` makes the haunch a plain frame: a VariableElement's
# stiffness is not Zygote-differentiable (it mutates), case loads or not
function case_model(; variable::Bool = true)
    shmat = Material(3.0e7, 2.4, 0.2)
    slab = ShellSection(shmat, 0.2)
    P = [0.0 2.0 4.0 0.0 2.0 4.0 0.0 2.0 4.0;
         0.0 0.0 0.0 2.0 2.0 2.0 4.0 4.0 4.0]
    cells = [(1, 2, 5, 4), (2, 3, 6, 5), (4, 5, 8), (4, 8, 7), (5, 6, 9), (5, 9, 8)]
    m = shell_mesh(P, cells, slab)
    n = m.nodes
    for k in (1, 3)
        fixnode!(n[k], Bool[0, 0, 0, 0, 0, 0])
    end
    base = Node([4.0, 4.0, -3.0], :fixed)
    fmat = Material(2.0e8, 7.7e7, 7.85, 0.3)
    sec = Section(fmat, 1e-2, 2e-4, 1e-4, 5e-5)
    big = Section(fmat, 2e-2, 6e-4, 2e-4, 1e-4)
    edge = FrameElement(n[7], n[8], sec, :edge)
    haunch = variable ? VariableElement(n[8], n[9], AbstractSection{Float64}[big, sec], [0.5], :haunch) :
             FrameElement(n[8], n[9], big, :haunch)
    col = FrameElement(base, n[9], sec, :col)
    quads = filter(el -> el isa ShellElement{Float64,4}, m.elements)
    tris = filter(el -> el isa ShellElement{Float64,3}, m.elements)
    loads = AbstractLoad{Float64}[
        [SelfWeight(el; g = [0.0, 0, -9.81], case = :dead) for el in m.elements]...,
        LineLoad(edge, [0.0, 0.0, -5.0]; case = :dead),
        LineLoad(haunch, [0.0, 0.0, -3.0]; case = :dead),          # VariableElement: static
        AreaLoad(quads[1], -2.0; case = :live), AreaLoad(tris[2], -3.0; case = :live),
        PointLoad(edge, 0.3, [0.0, 0.0, -10.0]; case = :live),
        NodeForce(n[5], [0.0, 0.0, -4.0]; case = :live),
        NodeForce(n[9], [6.0, 0.0, 0.0]; case = :wind),
        LineLoad(col, [1.5, 0.0, 0.0]; case = :wind)]
    return Model([n; base], AbstractElement{Float64}[m.elements..., edge, haunch, col], loads)
end

# the same state with every free node moved a little and every section scaled
function perturbed_state(model)
    s = extract_state(model)
    X = s.X .+ 0.03 .* sin.(reshape(1:length(s.X), size(s.X)))
    secs = map(s.sections) do sec
        sec isa ShellSection ? ShellSection(sec.material, 1.1 * sec.t) :
        sec isa Section ? Section(sec.material, 1.2 * sec.A, sec.Ix, sec.Iy, sec.J) : sec
    end
    return ModelState{Float64}(X, secs, nothing, nothing)
end

@testset "Load cases on the pure path (C0: CaseLoads)" begin
    model = case_model()
    process!(model)
    cache = model.cache

    @testset "table: names, sizes, show" begin
        cl = case_loads(model)
        @test cl.names == [:dead, :live, :wind]
        @test length(cl) == 3
        n = cache.partition.n_global
        @test size(cl.P) == size(cl.static) == size(cl.F0) == (n, 3)
        @test length(cl.fef.loads) == 3 && length(cl.fef.sloads) == 3 && length(cl.fef.tloads) == 5
        @test cl.fcase == [1, 2, 3]
        @test count(==(2), cl.scase) == 1 && count(==(2), cl.tcase) == 1
        @test size(cl.cscatter) == (3n, 36) && nnz(cl.cscatter) == nnz(cl.fef.scatter)
        @test occursin(":dead, :live, :wind", sprint(show, MIME"text/plain"(), cl))
        sub = case_loads(model; cases = [:wind, :dead])
        @test sub.names == [:wind, :dead] && sub.P[:, 2] == cl.P[:, 1]
        @test_throws ArgumentError case_loads(model; cases = [:snow])
        @test_throws ArgumentError case_loads(model; cases = [:dead, :dead])
        @test_throws ArgumentError case_loads(model; cases = Symbol[])
        @test_throws ArgumentError case_loads(case_model())                 # not processed
    end

    @testset "frozen load vectors ≡ solve_cases!" begin
        cl = case_loads(model)
        cr = solve_cases!(model)
        @test cr.cases == cl.names
        @test cl.F0 == cr.F
        # the static part is the VariableElement's line load, in :dead only
        @test norm(cl.static[:, 1]) > 0 && iszero(cl.static[:, 2:3])
    end

    @testset "the cache is left exactly as it was" begin
        assemble_loads!(cache, model; case = :live)          # an unusual prior state
        P, Pf, fef = copy(cache.P), copy(cache.Pf), cache.fef
        q = deepcopy(cache.q_local)
        case_loads(model)
        @test cache.P == P && cache.Pf == Pf && cache.fef === fef && cache.q_local == q
        assemble_loads!(cache, model)
    end

    @testset "design-dependent fixed-end forces per case" begin
        cl = case_loads(model)
        s0 = extract_state(model)
        # at the model's own state: the cached per-case fixed-end forces
        @test case_Pf(cl, s0) ≈ cl.P - cl.F0 rtol = 1e-12 atol = 1e-12
        # at another state: each column is the single-case assemble_Pf there
        s1 = perturbed_state(model)
        Pf1 = case_Pf(cl, s1)
        for (k, c) in enumerate(cl.names)
            assemble_loads!(cache, model; case = c)
            @test Pf1[:, k] ≈ assemble_Pf(cache, s1) rtol = 1e-12 atol = 1e-12
        end
        # summed over the cases: the all-loads fixed-end forces
        assemble_loads!(cache, model)
        @test vec(sum(Pf1; dims = 2)) ≈ assemble_Pf(cache, s1) rtol = 1e-12 atol = 1e-12
        @test norm(Pf1 - case_Pf(cl, s0)) > 1e-6                          # it did move
    end
end

# a rigid-floor storey (diaphragm) with two cases, one of them a column line load
function case_storey()
    mat = Material(2.0e8, 7.85, 0.3)
    xy = ((0.0, 0.0), (6.0, 0.0), (0.0, 4.0), (6.0, 4.0))
    base = [Node([x, y, 0.0], :fixed) for (x, y) in xy]
    top = [Node([x, y, 3.0], :free) for (x, y) in xy]
    cols = AbstractElement{Float64}[FrameElement(base[k], top[k], Section(mat, 1e-2, I, I, 0.5I), :col)
                                    for (k, I) in enumerate((1e-4, 2e-4, 1.5e-4, 3e-4))]
    d = RigidDiaphragm(top; id = :roof)
    loads = AbstractLoad{Float64}[NodeForce(d.master, [10.0, -4.0, 0.0]; case = :wind),
        NodeMoment(d.master, [0.0, 0.0, 25.0]; case = :wind),
        LineLoad(cols[2], [0.0, 2.0, 0.0]; case = :side)]
    return Model(vcat(base, top), cols, loads; diaphragms = [d])
end

@testset "Load cases on the pure path (C1: solve, compliance)" begin
    model = case_model()
    process!(model)
    cache = model.cache
    cl = case_loads(model)
    s0, s1 = extract_state(model), perturbed_state(model)

    @testset "at the model's state ≡ solve_cases!" begin
        cr = solve_cases!(model)
        U = solve(model, s0, cl)
        @test size(U) == (cache.partition.n_global, 3)
        for k in 1:3
            @test U[:, k] ≈ cr.results[k].u rtol = 1e-10
        end
        @test compliance(model, s0, cl) ≈ [r.compliance for r in cr.results] rtol = 1e-10
        @test solve(model, s0, cl; loads = :frozen) ≈ U rtol = 1e-10
        # combinations: U·λ and (Fλ)ᵀ(Uλ) are combine()'s displacements and compliance
        combo = LoadCombination(:LRFD, [:dead => 1.2, :live => 1.6, :wind => 0.5])
        λ = combination_weights(cl, combo)
        @test λ == [1.2, 1.6, 0.5]
        rc = combine(cr, combo)
        @test U * λ ≈ rc.u rtol = 1e-10
        free = cache.partition.free
        @test dot((cl.P - case_Pf(cl, s0))[free, :] * λ, (U * λ)[free]) ≈ rc.compliance rtol = 1e-10
        @test combination_weights(cl, LoadCombination(:svc, [:dead => 1.0])) == [1.0, 0.0, 0.0]
        @test_throws ArgumentError combination_weights(cl, LoadCombination(:x, [:snow => 1.0]))
    end

    @testset "at another state ≡ one single-case solve per case" begin
        U1 = solve(model, s1, cl)
        U1f = solve(model, s1, cl; loads = :frozen)
        c1 = compliance(model, s1, cl)
        for (k, c) in enumerate(cl.names)
            assemble_loads!(cache, model; case = c)
            @test U1[:, k] ≈ solve(model, s1) rtol = 1e-10
            @test U1f[:, k] ≈ solve(model, s1; loads = :frozen) rtol = 1e-10
            @test c1[k] ≈ compliance(model, s1) rtol = 1e-10
        end
        assemble_loads!(cache, model)
        @test vec(sum(U1; dims = 2)) ≈ solve(model, s1) rtol = 1e-10        # superposition
        @test norm(U1 - U1f) > 1e-8 * norm(U1)                             # design loads moved
        @test_throws ArgumentError solve(model, s1, cl; loads = :fixed)
    end

    @testset "rigid diaphragm" begin
        m = case_storey()
        solve!(m)
        c = case_loads(m)
        s = extract_state(m)
        s2 = ModelState{Float64}(s.X .+ 0.02 .* cos.(reshape(1:length(s.X), size(s.X))), s.sections, nothing, nothing)
        U = solve(m, s2, c)
        for (k, name) in enumerate(c.names)
            assemble_loads!(m.cache, m; case = name)
            @test U[:, k] ≈ solve(m, s2) rtol = 1e-10
        end
        assemble_loads!(m.cache, m)
    end

    @testset "gradients: Zygote ≡ ForwardDiff ≡ finite differences" begin
        model = case_model(; variable = false)
        process!(model)
        cl = case_loads(model)
        s0, s1 = extract_state(model), perturbed_state(model)
        cfg = Zygote.ZygoteRuleConfig()
        @test rrule(cfg, Asap._frame_fef_q, cl.fef, s0) !== nothing
        @test rrule(cfg, Asap._shell_fef_q, cl.fef, s0) !== nothing
        @test rrule(cfg, Asap._shell3_fef_q, cl.fef, s0) !== nothing
        nX = length(s0.X)
        shells = findall(sec -> sec isa ShellSection, s0.sections)
        frames = findall(sec -> sec isa Section, s0.sections)
        # design: every node position, every shell thickness, every frame area
        function state(x)
            X = reshape(x[1:nX], 3, :)
            secs = map(eachindex(s0.sections)) do e
                sec = s0.sections[e]
                i = findfirst(==(e), shells)
                i !== nothing && return ShellSection(sec.material, x[nX+i])
                j = findfirst(==(e), frames)
                j !== nothing && return Section(sec.material, x[nX+length(shells)+j], sec.Ix, sec.Iy, sec.J)
                return sec
            end
            return ModelState{eltype(x)}(X, secs, nothing, nothing)
        end
        λ = [1.2, 1.6, 0.5]
        f(x) = dot(λ, compliance(model, state(x), cl))
        x0 = vcat(vec(s1.X), [s1.sections[e].t for e in shells], [s1.sections[e].A for e in frames])
        gz = Zygote.gradient(f, x0)[1]
        gf = ForwardDiff.gradient(f, x0)
        # finite differences along one direction (cheap; every entry enters)
        d = sin.(1:length(x0)) .* max.(abs.(x0), 1e-2)
        dd = central_fdm(5, 1)(ε -> f(x0 + ε * d), 0.0)
        @test gz ≈ gf rtol = 1e-8
        @test dot(gz, d) ≈ dd rtol = 1e-6
        @test norm(gz[nX+1:end]) > 0
    end
end

@testset "Load cases on the pure path (C2: per-case local fixed-end forces)" begin
    model = case_model(; variable = false)
    process!(model)
    cache = model.cache
    cl = case_loads(model)
    s1 = perturbed_state(model)
    Q = fef_local(cl, s1)
    @test size(Q) == (12, length(cl.fef.loads))
    # each case's columns ≡ the single-case fef_local of that case at the state
    for (k, c) in enumerate(cl.names)
        assemble_loads!(cache, model; case = c)
        @test Q[:, cl.fcase .== k] ≈ fef_local(cache, s1) rtol = 1e-12
        @test cl.fef.elem[cl.fcase .== k] == cache.fef.elem
    end
    assemble_loads!(cache, model)
    @test size(fef_local(case_loads(model; cases = [:live]), s1), 1) == 12
    # differentiable: Zygote ≡ ForwardDiff for a scalar of the per-case axial terms
    nX = length(s1.X)
    w = collect(1.0:size(Q, 2))
    f(x) = dot(w, fef_local(cl, ModelState{eltype(x)}(reshape(x, 3, :), s1.sections, nothing, nothing))[7, :])
    x0 = vec(s1.X)
    @test Zygote.gradient(f, x0)[1] ≈ ForwardDiff.gradient(f, x0) rtol = 1e-8
    @test rrule(Zygote.ZygoteRuleConfig(), fef_local, cl, s1) !== nothing
end
