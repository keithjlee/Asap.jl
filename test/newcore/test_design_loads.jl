# Design-dependent element loads on the pure path (v1.2): the fixed-end
# forces are re-evaluated at the state (assemble_Pf), so
#   - VALUE: the pure solve at any state equals a rebuilt model's stateful
#     solve (the ground truth), for every element-load type, under changes of
#     geometry, sections (ρA, EI) and end springs;
#   - GRADIENT: Zygote and ForwardDiff through the pure path agree with
#     central differences of the REBUILT model's compliance — the check that
#     caught frozen loads (AsapOptim gradients were off by −35× … −0.008×).

using ForwardDiff

const DL_STEEL = AsapNext.Material(200e6, 80e6, 8.0, 0.3)

# a 3D three-member frame carrying every element-load type: line load,
# point load, point moment, self-weight, partial-span local distributed load,
# with rolled sections, a released end and semi-rigid beam ends
function dl_frame(X; A = 4e-3, Ix = 2e-5, k = 1e4)
    N = AsapNext
    s = N.Section(DL_STEEL, A, Ix, 6e-6, 1e-7)
    n = [N.Node(X[:, 1], :fixed), N.Node(X[:, 2], :free), N.Node(X[:, 3], :free), N.Node(X[:, 4], :fixed)]
    sp = N.EndSprings(Inf, Inf, k, k)
    els = N.AbstractElement{Float64}[N.FrameElement(n[1], n[2], s; rollangle = 0.2),
        N.FrameElement(n[2], n[3], s, N.EndConditions(sp, sp), :beam),
        N.FrameElement(n[3], n[4], s; release = :fixedfree, rollangle = 1.1)]
    loads = N.AbstractLoad{Float64}[N.LineLoad(els[2], [0.0, -5.0, 1.0]),
        N.PointLoad(els[1], 0.4, [1.0, 0.0, -2.0]),
        N.PointMoment(els[3], 0.3, [0.5, 1.0, -1.0]), N.SelfWeight(els[2]),
        N.DistributedLoad(els[3], [0.1, 0.6], [2.0, -1.0], [0.0, 0.0, 1.0]; coords = :local),
        N.NodeForce(n[2], [3.0, 0.0, 0.0])]
    return N.Model(n, els, loads)
end

const DL_X0 = [0.0 0.0 4.0 4.5; 0.0 3.0 3.2 0.0; 0.0 0.5 0.2 0.0]
const DL_DX = [0.0 0.3 -0.2 0.0; 0.0 0.4 0.1 0.0; 0.0 -0.2 0.3 0.0]

# compliance of a rebuilt, stateful-solved model (the ground truth)
function dl_true_compliance(model)
    N = AsapNext
    N.solve!(model)
    F = model.cache.P - model.cache.Pf
    return sum(F .* model.results.u)
end

@testset "design-dependent element loads" begin
    N = AsapNext
    m0 = dl_frame(DL_X0)
    N.process!(m0)
    N.assemble_loads!(m0.cache, m0)
    st0 = N.extract_state(m0)

    @testset "pure kernels reproduce the stateful assembly" begin
        @test N.assemble_Pf(m0.cache, st0) ≈ m0.cache.Pf rtol = 1e-14
        Q = N.fef_local(m0.cache, st0)
        @test size(Q) == (12, 5)
        # per-element sum of the per-load q̃ = the stateful q_local
        for e in eachindex(m0.elements)
            cols = findall(==(e), m0.cache.fef.elem)
            @test vec(sum(Q[:, cols]; dims = 2)) ≈ m0.cache.q_local[e][1] atol = 1e-12
        end
        # :frozen is the cached Pf, exactly
        @test N.solve(m0, st0; loads = :frozen) ≈ N.solve(m0, st0) rtol = 1e-12
        @test_throws ArgumentError N.solve(m0, st0; loads = :bogus)
    end

    @testset "value parity with rebuilt models" begin
        # geometry: every load type rides with its member
        X1 = DL_X0 + DL_DX
        m1 = dl_frame(X1)
        N.solve!(m1)
        st1 = N.ModelState{Float64}(X1, st0.sections, st0.EA, nothing)
        @test N.solve(m0, st1) ≈ m1.results.u rtol = 1e-12
        # … and the frozen path is genuinely different (the old behavior)
        @test norm(N.solve(m0, st1; loads = :frozen) - m1.results.u) > 1e-3 * norm(m1.results.u)

        # sections (ρA → self-weight, EI → semi-rigid condensation) and end springs
        m2 = dl_frame(DL_X0; A = 6e-3, Ix = 3e-5, k = 3e3)
        N.solve!(m2)
        s2 = N.extract_state(m2)
        st2 = N.ModelState{Float64}(DL_X0, s2.sections, s2.EA, Any[el.ends for el in m2.elements])
        @test N.solve(m0, st2) ≈ m2.results.u rtol = 1e-12
        @test N.compliance(m0, st2) ≈ dl_true_compliance(m2) rtol = 1e-12
    end

    @testset "gradients: geometry (Zygote, ForwardDiff vs rebuilt FD)" begin
        obj(x) = N.compliance(m0, N.ModelState{eltype(x)}(reshape(x, 3, :), st0.sections))
        x0 = vec(DL_X0)
        free = [4, 5, 6, 7, 8, 9]                      # coordinates of nodes 2, 3
        g_zy = Zygote.gradient(obj, x0)[1]
        g_fd = ForwardDiff.gradient(obj, x0)
        g_true = map(free) do i
            h = 1e-6
            Xp = copy(DL_X0); Xp[i] += h
            Xm = copy(DL_X0); Xm[i] -= h
            (dl_true_compliance(dl_frame(Xp)) - dl_true_compliance(dl_frame(Xm))) / 2h
        end
        @test g_zy[free] ≈ g_true rtol = 1e-5
        @test g_fd[free] ≈ g_true rtol = 1e-5
        @test g_zy ≈ g_fd rtol = 1e-10
    end

    @testset "gradients: area via self-weight, inertia via condensation" begin
        function obj(p)       # p = (A, Ix) of every element
            secs = map(i -> N.Section(DL_STEEL, p[1], p[2], 6e-6, 1e-7), 1:3)
            return N.compliance(m0, N.ModelState{eltype(p)}(DL_X0, secs))
        end
        p0 = [4e-3, 2e-5]
        g_zy = Zygote.gradient(obj, p0)[1]
        g_fd = ForwardDiff.gradient(obj, p0)
        hA, hI = 1e-8, 1e-8           # (smaller steps drown in solver roundoff: C ≈ 10)
        g_true = [(dl_true_compliance(dl_frame(DL_X0; A = p0[1] + hA)) -
                   dl_true_compliance(dl_frame(DL_X0; A = p0[1] - hA))) / 2hA,
                  (dl_true_compliance(dl_frame(DL_X0; Ix = p0[2] + hI)) -
                   dl_true_compliance(dl_frame(DL_X0; Ix = p0[2] - hI))) / 2hI]
        @test g_zy ≈ g_true rtol = 1e-4
        @test g_fd ≈ g_true rtol = 1e-4
        @test g_zy ≈ g_fd rtol = 1e-10
    end

    @testset "gradients: joint stiffness with a loaded beam" begin
        # the case AsapOptim used to forbid: κ moves the condensed end moments
        # of the loaded beam, from 0 (hinged) toward wL²/12 (clamped)
        function obj(k)
            sp = N.EndSprings(Inf, Inf, k[1], k[1])
            ends = Any[m0.elements[1].ends, N.EndConditions(sp, sp), m0.elements[3].ends]
            return N.compliance(m0, N.ModelState{eltype(k)}(DL_X0, st0.sections, nothing, ends))
        end
        k0 = [1e4]
        g_zy = Zygote.gradient(obj, k0)[1]
        g_fd = ForwardDiff.gradient(obj, k0)
        h = 1.0
        g_true = (dl_true_compliance(dl_frame(DL_X0; k = k0[1] + h)) -
                  dl_true_compliance(dl_frame(DL_X0; k = k0[1] - h))) / 2h
        @test g_zy[1] ≈ g_true rtol = 1e-5
        @test g_fd[1] ≈ g_true rtol = 1e-5
    end

    @testset "kernel semantics under design change" begin
        # a uniform line load's total scales with the CURRENT length
        mat = DL_STEEL
        sec = N.Section(mat, 4e-3, 2e-5, 6e-6, 1e-7)
        el = N.FrameElement(N.Node([0.0, 0, 0], :free), N.Node([1.0, 0, 0], :free), sec)
        w = 3.0
        load = N.LineLoad(el, [0.0, -w, 0.0])
        for L in (1.0, 2.5, 7.0)
            q = N.fixed_end_forces(load, sec, el.ends, [0.0, 0, 0], [L, 0, 0], el.rollangle)
            @test abs(sum(q[[2, 3, 8, 9]])) ≈ w * L rtol = 1e-12   # |total transverse FEF| = wL
        end
        # a point load stays at its fraction of the member
        pl = N.PointLoad(el, 0.25, [0.0, -1.0, 0.0])
        for L in (1.0, 4.0)
            q = N.fixed_end_forces(pl, sec, el.ends, [0.0, 0, 0], [L, 0, 0], el.rollangle)
            # clamped end moment magnitudes P·a·b²/L² at a = L/4
            a, b = L / 4, 3L / 4
            @test maximum(abs, q[[5, 6]]) ≈ a * b^2 / L^2 rtol = 1e-12
        end
        # condensation promotes: a Dual end spring (or EI) with a plain-Float64
        # clamped q — the implicit-diff seed pattern — must not truncate
        kd = ForwardDiff.Dual(1e4, 1.0)
        spd = N.EndSprings(Inf, Inf, kd, kd)
        q = N.fixed_end_forces(load, sec, el.ends, [0.0, 0, 0], [3.0, 0, 0], el.rollangle)
        q̃ = N.condense_fef(q, sec, 3.0, N.EndConditions(spd, spd))
        @test eltype(q̃) <: ForwardDiff.Dual
        @test any(!iszero, ForwardDiff.partials.(q̃, 1))
        # self-weight differentiates in ρA (ForwardDiff Dual area)
        sw = N.SelfWeight(el)
        f(A) = sum(N.fixed_end_forces(sw, N.Section(mat, A, 2e-5, 6e-6, 1e-7), el.ends,
            [0.0, 0, 0], [2.0, 0, 0], el.rollangle))
        @test ForwardDiff.derivative(f, 4e-3) ≈ (f(4e-3 + 1e-9) - f(4e-3 - 1e-9)) / 2e-9 rtol = 1e-6
    end
end
