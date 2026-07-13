# Phase 4: load cases, combinations, envelopes.
# Exit criteria: combination results ≡ brute-force solve of the factored
# loads; envelope ≡ station-wise extrema of brute-force per-combo diagrams;
# all from a single factorization.

# a portal frame with three tagged cases: dead (self-weight-ish line load),
# live (point + line), wind (lateral)
function case_portal(; factored::Union{Nothing,Dict{Symbol,Float64}}=nothing)
    N = AsapNext
    mat = N.Material(200.0, 77.0, 1.0, 0.3)
    sec = N.Section(mat, 1e4, 8e7, 3e7, 5e6)
    H, L = 3000.0, 5000.0
    n1 = N.Node([0.0, 0.0, 0.0], :fixed)
    n2 = N.Node([0.0, 0.0, H], :free)
    n3 = N.Node([L, 0.0, H], :free)
    n4 = N.Node([L, 0.0, 0.0], :fixed)
    beam = N.FrameElement(n2, n3, sec, :beam)
    els = N.AbstractElement{Float64}[
        N.FrameElement(n1, n2, sec, :col), beam, N.FrameElement(n4, n3, sec, :col)]

    λ(c) = factored === nothing ? 1.0 : get(factored, c, 0.0)
    loads = N.AbstractLoad{Float64}[]
    λ(:dead) != 0 && push!(loads,
        N.LineLoad(beam, λ(:dead) * [0.0, 0.0, -2.0]; case=:dead))
    λ(:live) != 0 && append!(loads, [
        N.PointLoad(beam, 0.4, λ(:live) * [0.0, 0.0, -10.0]; case=:live),
        N.LineLoad(beam, λ(:live) * [0.0, 0.0, -1.0]; case=:live)])
    λ(:wind) != 0 && push!(loads,
        N.NodeForce(n2, λ(:wind) * [5.0, 0.0, 0.0]; case=:wind))

    return N.Model([n1, n2, n3, n4], els, loads)
end

@testset "Cases, combinations, envelopes (Phase 4)" begin
    N = AsapNext

    model = case_portal()
    cr = N.solve_cases!(model)

    @testset "case bookkeeping" begin
        @test cr.cases == [:dead, :live, :wind]
        @test N.load_cases(model) == [:dead, :live, :wind]
        # each per-case result reflects only its own loads: wind causes no
        # vertical beam deflection at midspan-symmetric points, dead no sway
        wind = cr[:wind]
        dead = cr[:dead]
        n2, n3 = model.nodes[2], model.nodes[3]
        # symmetric gravity spreads the frame: top nodes move apart
        # (antisymmetric u_x); wind sways it: both move the same way
        @test N.displacement(dead, n2)[1] ≈ -N.displacement(dead, n3)[1] rtol = 1e-9
        @test sign(N.displacement(wind, n2)[1]) == sign(N.displacement(wind, n3)[1])
        @test abs(N.displacement(wind, n2)[1]) > 0
        # gravity dominates vertical response
        @test abs(N.displacement(dead, n2)[3]) > abs(N.displacement(wind, n2)[3])
    end

    combos = [
        N.LoadCombination(:strength, [:dead => 1.2, :live => 1.6, :wind => 0.5]),
        N.LoadCombination(:service, [:dead => 1.0, :live => 1.0]),
        N.LoadCombination(:uplift, [:dead => 0.9, :wind => 1.6]),
    ]

    @testset "combine ≡ brute-force factored solve: $(c.name)" for c in combos
        res = N.combine(cr, c)

        brute_model = case_portal(; factored=N.factor_dict(c))
        N.solve!(brute_model)
        brute = brute_model.results

        @test res.u ≈ brute.u rtol = 1e-10 atol = 1e-12
        @test res.reactions ≈ brute.reactions rtol = 1e-10 atol = 1e-10
        for i in eachindex(res.element_forces)
            @test res.element_forces[i] ≈ brute.element_forces[i] rtol = 1e-9 atol = 1e-8
        end
        @test res.compliance ≈ brute.compliance rtol = 1e-9

        # combination-aware recovery matches the brute model's diagrams
        beam = model.elements[2]
        bbeam = brute_model.elements[2]
        st = N.internal_forces(model, beam; results=res, factors=N.factor_dict(c))
        stb = N.internal_forces(brute_model, bbeam)
        for t in (0.1, 0.35, 0.62, 0.9)
            @test N.moment_z(st, t) ≈ N.moment_z(stb, t) rtol = 1e-8 atol = 1e-6
            @test N.moment_y(st, t) ≈ N.moment_y(stb, t) rtol = 1e-8 atol = 1e-6
            @test N.shear_z(st, t) ≈ N.shear_z(stb, t) rtol = 1e-8 atol = 1e-8
            @test N.axial_force(st, t) ≈ N.axial_force(stb, t) rtol = 1e-8 atol = 1e-8
        end
    end

    @testset "envelope ≡ brute-force extrema" begin
        beam = model.elements[2]
        env = N.envelope(model, beam, cr, combos; resolution=15)
        @test env.combos == [:strength, :service, :uplift]
        @test all(env.lo .<= env.hi)

        # brute force: per-combo factored model diagrams at the same stations
        for (q, row) in ((:My, 5), (:Vz, 4), (:N, 1))
            per_combo = map(combos) do c
                bm = case_portal(; factored=N.factor_dict(c))
                N.solve!(bm)
                stb = N.internal_forces(bm, bm.elements[2])
                fn = Dict(:N => N.axial_force, :Vz => N.shear_z, :My => N.moment_y)[q]
                [fn(stb, xi / stb.L) for xi in env.x]
            end
            brute_lo = [minimum(pc[k] for pc in per_combo) for k in eachindex(per_combo[1])]
            brute_hi = [maximum(pc[k] for pc in per_combo) for k in eachindex(per_combo[1])]
            @test env.lo[row, :] ≈ brute_lo rtol = 1e-8 atol = 1e-6
            @test env.hi[row, :] ≈ brute_hi rtol = 1e-8 atol = 1e-6
        end
    end

    @testset "single factorization reused" begin
        # after solve_cases!, the cache factorization solves any case RHS to
        # the already-computed displacements (proof the factorization spans
        # all cases; a fresh factorize would be a different object)
        cache = model.cache
        fact = cache.factorization
        free = cache.partition.free
        for (j, case) in enumerate(cr.cases)
            uf = fact \ cr.F[free, j]
            @test uf ≈ cr.results[j].u[free] rtol = 1e-12
        end
    end
end
