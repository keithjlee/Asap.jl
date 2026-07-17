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
