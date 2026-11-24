# Shell S6: the differentiable pure path with shells — parity with the
# in-place solve, design-dependent shell loads (area loads follow the
# element area/orientation, self-weight follows ρt), and gradients w.r.t.
# node positions and thickness against finite differences (Zygote). The fast
# reverse rules live in ext/AsapShellADExt.jl (needs ChainRulesCore + ForwardDiff).

using StaticArrays, LinearAlgebra, Zygote, FiniteDifferences
import ForwardDiff

function _ad_wall(; with_diaphragm = false)
    mat = Material(3.0e4, 2.4e-3, 0.2)
    sec = ShellSection(mat, 0.2)
    w = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 3.0, 2.0, sec;
        openings = [(1.0, 0.0, 0.8, 1.2)], size = 0.5)
    loads = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81]) for el in w.elements]...,
                                  [AreaLoad(el, 0.05) for el in w.elements]...]    # out-of-plane pressure
    top = [w.grid[i, end] for i in eachindex(w.s)]
    append!(loads, [NodeForce(n, [2.0 / length(top), 0, 0]) for n in top])
    ds = with_diaphragm ? [RigidDiaphragm(top)] : RigidDiaphragm{Float64}[]
    m = Model(w.nodes, AbstractElement{Float64}[w.elements...], loads; diaphragms = ds)
    solve!(m)
    return m, w, mat
end

@testset "Shell AD / pure path (S6)" begin

    @testset "pure path ≡ in-place (design-dependent shell loads)" begin
        for dia in (false, true)
            m, w, mat = _ad_wall(; with_diaphragm = dia)
            st = extract_state(m)
            @test compliance(m, st) ≈ m.results.compliance rtol = 1e-10
            @test compliance(m, st; loads = :frozen) ≈ m.results.compliance rtol = 1e-10
            @test solve(m, st) ≈ m.results.u rtol = 1e-9
            @test Asap.assemble_Pf(m.cache, st) ≈ m.cache.Pf rtol = 1e-12
        end
    end

    @testset "gradients w.r.t. node positions" begin
        m, w, mat = _ad_wall()
        st = extract_state(m)
        f(X) = compliance(m, ModelState{eltype(X)}(X, st.sections))
        gz = Zygote.gradient(f, st.X)[1]
        gfd = FiniteDifferences.grad(central_fdm(5, 1), f, st.X)[1]
        @test gz ≈ gfd rtol = 1e-5
        # (full-pipeline ForwardDiff ≡ Zygote to 1e-8 was verified in development;
        # left out of the suite — its first compile through the shell kernels
        # costs ~100 s)
    end

    @testset "gradient w.r.t. thickness (stiffness AND self-weight)" begin
        m, w, mat = _ad_wall()
        st = extract_state(m)
        idx = [findfirst(e -> e === el, m.elements) for el in w.elements]
        function f(t)
            secs = Any[st.sections...]
            secs = map(eachindex(secs)) do i
                i in idx ? ShellSection(mat, t[1]) : secs[i]
            end
            compliance(m, ModelState{Float64}(st.X, secs))
        end
        t0 = [0.2]
        gz = Zygote.gradient(f, t0)[1]
        gfd = FiniteDifferences.grad(central_fdm(5, 1), f, t0)[1]
        @test gz ≈ gfd rtol = 1e-6
        # frozen loads drop the self-weight term: a different (wrong for
        # thickness design) gradient — the design-dependent path matters
        ffz(t) = compliance(m, ModelState{Float64}(st.X,
            map(i -> i in idx ? ShellSection(mat, t[1]) : st.sections[i], eachindex(st.sections))); loads = :frozen)
        @test !(Zygote.gradient(ffz, t0)[1] ≈ gz)
    end
end
