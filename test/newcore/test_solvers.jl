#=
Solver seam + LinearSolve.jl backends: every backend must reproduce the
built-in CHOLMOD path, stay sticky on the cache, and survive
mutate-and-resolve refactorization.
=#
using LinearSolve

@testset "solver backends (LinearSolve extension)" begin
    N = AsapNext
    sec = N.Section(N.Material(200e6, 77e6, 8.0, 0.3), 1e-2, 1e-4, 5e-5, 1e-6)

    function portal()
        n1 = N.Node([0.0, 0.0, 0.0], :fixed)
        n2 = N.Node([0.0, 0.0, 3.0], :free)
        n3 = N.Node([5.0, 0.0, 3.0], :free)
        n4 = N.Node([5.0, 0.0, 0.0], :fixed)
        b = N.FrameElement(n2, n3, sec, :beam)
        els = N.AbstractElement{Float64}[N.FrameElement(n1, n2, sec), b, N.FrameElement(n4, n3, sec)]
        loads = N.AbstractLoad{Float64}[N.LineLoad(b, [0.0, 0.0, -2.0]),
            N.NodeForce(n2, [5.0, 0.0, 0.0])]
        N.Model([n1, n2, n3, n4], els, loads)
    end

    ref = portal()
    N.solve!(ref)                                  # built-in CHOLMOD path
    uref = copy(ref.results.u)

    @testset "$(nameof(typeof(alg)))" for alg in (KLUFactorization(), KrylovJL_CG())
        m = portal()
        N.solve!(m; solver = alg)
        @test m.cache.factorization isa N.FactorizationCache
        @test m.cache.factorization.solver === alg
        @test m.results.u ≈ uref rtol = 1e-8

        # sticky: subsequent solves keep the backend, refactorize in place
        fc = m.cache.factorization
        m.elements[1].section = N.Section(N.Material(200e6, 77e6, 8.0, 0.3), 2e-2, 2e-4, 1e-4, 2e-6)
        N.solve!(m)
        @test m.cache.factorization === fc

        ref2 = portal()
        ref2.elements[1].section = m.elements[1].section
        N.solve!(ref2)
        @test m.results.u ≈ ref2.results.u rtol = 1e-8
    end

    @testset "solve_cases! through a backend" begin
        m = portal()
        m.loads[1] = N.LineLoad(m.elements[2], [0.0, 0.0, -2.0]; case = :dead)
        m.loads[2] = N.NodeForce(m.nodes[2], [5.0, 0.0, 0.0]; case = :wind)
        cr = N.solve_cases!(m; solver = KLUFactorization())
        mref = portal()
        mref.loads[1] = N.LineLoad(mref.elements[2], [0.0, 0.0, -2.0]; case = :dead)
        mref.loads[2] = N.NodeForce(mref.nodes[2], [5.0, 0.0, 0.0]; case = :wind)
        crref = N.solve_cases!(mref)
        for c in (:dead, :wind)
            @test cr[c].u ≈ crref[c].u rtol = 1e-8
        end
    end

    @testset "switching backends rebuilds the factorization" begin
        m = portal()
        N.solve!(m)
        fc0 = m.cache.factorization
        N.solve!(m; solver = KLUFactorization())
        @test m.cache.factorization !== fc0
        @test m.results.u ≈ uref rtol = 1e-8
    end
end

@testset "solver seam on the PURE path (solve_free 3-arg)" begin
    N = AsapNext
    sec = N.Section(N.Material(200e6, 77e6, 8.0, 0.3), 1e-2, 1e-4, 5e-5, 1e-6)
    n1 = N.Node([0.0, 0.0, 0.0], :fixed)
    n2 = N.Node([0.0, 0.0, 3.0], :free)
    n3 = N.Node([5.0, 0.0, 3.0], :free)
    n4 = N.Node([5.0, 0.0, 0.0], :fixed)
    els = N.AbstractElement{Float64}[N.FrameElement(n1, n2, sec),
        N.FrameElement(n2, n3, sec), N.FrameElement(n4, n3, sec)]
    model = N.Model([n1, n2, n3, n4], els,
        N.AbstractLoad{Float64}[N.NodeForce(n2, [5.0, 0.0, -3.0])])
    N.solve!(model)
    st = N.extract_state(model)

    u0 = N.solve(model, st)                                     # default backend
    for solver in (KLUFactorization(), KrylovJL_CG(), N.CachedSolver())
        @test N.solve(model, st; solver = solver) ≈ u0 rtol = 1e-8
    end

    # gradients flow identically through every backend (the 3-arg rrule)
    δ = zeros(size(st.X)); δ[3, 2] = 1.0
    obj(θ, solver) = N.compliance(model,
        N.ModelState{eltype(θ)}(st.X + θ[1] * δ, st.sections, st.EA, nothing);
        solver = solver)
    g0 = Zygote.gradient(θ -> obj(θ, nothing), [0.1])[1]
    for solver in (KLUFactorization(), N.CachedSolver())
        @test Zygote.gradient(θ -> obj(θ, solver), [0.1])[1] ≈ g0 rtol = 1e-8
    end
    gfd = FiniteDifferences.grad(central_fdm(5, 1), θ -> obj(θ, nothing), [0.1])[1]
    @test g0 ≈ gfd rtol = 1e-6

    # ForwardDiff through the 3-arg Dual methods
    for solver in (KLUFactorization(), N.CachedSolver())
        @test ForwardDiff.gradient(θ -> obj(θ, solver), [0.1]) ≈ g0 rtol = 1e-8
    end
end

@testset "CachedSolver reuse semantics" begin
    N = AsapNext
    cs = N.CachedSolver()
    # triplets in (1,1),(2,1),(2,2),(1,2) order — genuinely symmetric SPD
    K1 = sparse([1, 2, 2, 1], [1, 1, 2, 2], [4.0, 1.0, 3.0, 1.0])
    F = [1.0, 2.0]

    u1 = N.solve_free(cs, K1, F)
    @test u1 ≈ K1 \ F rtol = 1e-12
    @test cs.hits == 0 && cs.refactorizations == 0

    u1b = N.solve_free(cs, K1, F)                 # identical values → cache hit
    @test u1b ≈ u1
    @test cs.hits == 1 && cs.refactorizations == 0

    K2 = sparse([1, 2, 2, 1], [1, 1, 2, 2], [5.0, 0.5, 4.0, 0.5])  # same pattern
    u2 = N.solve_free(cs, K2, F)                  # numeric-only refactorization
    @test u2 ≈ K2 \ F rtol = 1e-12
    @test cs.refactorizations == 1

    u2b = N.solve_free(cs, K2, F)
    @test u2b ≈ u2 && cs.hits == 2
end
