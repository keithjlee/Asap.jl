# Phase 2 load types: trapezoids/triangles against classic closed-form FEF
# tables and an independent high-resolution quadrature; PointMoment against
# textbook concentrated-moment fixed-end values; case tags.

# independent oracle: EXACT closed-form integration of −N(x)ᵀw(x) for a
# linear w over [a, b] — polynomial antiderivatives, no quadrature, no code
# shared with the Gauss engine.
polyint(c, a, b) = sum(c[k] * (b^k - a^k) / k for k in eachindex(c))

function exact_trapezoid_fef(t1, t2, w1, w2, d, L)
    a, b = t1 * L, t2 * L
    # w(x) = α + βx on [a, b]
    β = (w2 - w1) / (b - a)
    α = w1 - β * a
    # shape functions as polynomial coefficients in x (constant term first)
    Na = [1.0, -1 / L]
    Nb = [0.0, 1 / L]
    H1 = [1.0, 0.0, -3 / L^2, 2 / L^3]
    H2 = [0.0, 1.0, -2 / L, 1 / L^2]
    H3 = [0.0, 0.0, 3 / L^2, -2 / L^3]
    H4 = [0.0, 0.0, -1 / L, 1 / L^2]
    # ∫ (α + βx)·p(x) dx = α·∫p + β·∫x·p  (x·p shifts coefficients up one power)
    I(p) = α * polyint(p, a, b) + β * polyint(vcat(0.0, p), a, b)

    q = zeros(12)
    q[1] = -d[1] * I(Na); q[7] = -d[1] * I(Nb)
    q[2] = -d[2] * I(H1); q[6] = -d[2] * I(H2)
    q[8] = -d[2] * I(H3); q[12] = -d[2] * I(H4)
    q[3] = -d[3] * I(H1); q[5] = d[3] * I(H2)
    q[9] = -d[3] * I(H3); q[11] = d[3] * I(H4)
    return q
end

@testset "Phase 2 loads" begin
    N = AsapNext
    mat = N.Material(200.0, 77.0, 1.0, 0.3)
    sec = N.Section(mat, 1e4, 8e7, 3e7, 5e6)
    L = 4000.0
    n1 = N.Node([0.0, 0.0, 0.0], :fixed)
    n2 = N.Node([L, 0.0, 0.0], :fixed)
    el = N.FrameElement(n1, n2, sec; rollangle=0.0)   # local y = global Y
    ends = N.EndConditions(:fixedfixed)
    x1, x2 = n1.position, n2.position

    fef(load) = N.fixed_end_forces(load, sec, ends, x1, x2, 0.0)

    @testset "full-span triangle ≡ classic table" begin
        # downward ramp 0 → W at the far end: |M1| = WL²/30, |M2| = WL²/20,
        # |V1| = 3WL/20, |V2| = 7WL/20 (Roark / AISC fixed-end tables)
        W = 5.0
        q = fef(N.TrapezoidLoad(el, 0.0, 1.0, 0.0, W, [0.0, -1.0, 0.0]))
        @test q[2] ≈ 3W * L / 20 rtol = 1e-12      # engine sign: reaction sense
        @test q[8] ≈ 7W * L / 20 rtol = 1e-12
        @test q[6] ≈ W * L^2 / 30 rtol = 1e-12
        @test q[12] ≈ -W * L^2 / 20 rtol = 1e-12
    end

    @testset "partial trapezoids ≡ independent quadrature" begin
        cases = [
            (0.0, 1.0, 2.0, 5.0),      # full-span trapezoid
            (0.25, 0.8, 0.0, 4.0),     # partial triangle
            (0.3, 0.6, 3.0, 1.0),      # partial decreasing trapezoid
        ]
        d = normalize([0.3, -0.8, 0.5])
        for (t1, t2, w1, w2) in cases
            q = fef(N.TrapezoidLoad(el, t1, t2, w1, w2, d))
            # direction is global; rollangle=0 element along X ⇒ local ≡ global here
            @test Vector(q) ≈ exact_trapezoid_fef(t1, t2, w1, w2, d, L) rtol = 1e-11
        end
    end

    @testset "PointMoment ≡ classic table" begin
        # concentrated moment M0 about local z at midspan of a clamped beam:
        # fixed-end moments M0/4 at both ends, shears 6·M0·ab/L³ = 1.5·M0/L
        M0 = 1000.0
        q = fef(N.PointMoment(el, 0.5, [0.0, 0.0, M0]))
        @test q[6] ≈ M0 / 4 rtol = 1e-12
        @test q[12] ≈ M0 / 4 rtol = 1e-12
        # shear couple magnitude 6·M0·ab/L³ = 1.5·M0/L (sign per the engine's
        # q = −consistent-load convention, verified by exact equilibrium below)
        @test q[2] ≈ 1.5 * M0 / L rtol = 1e-12
        @test q[8] ≈ -1.5 * M0 / L rtol = 1e-12
        @test q[2] + q[8] ≈ 0 atol = 1e-12
        # moment equilibrium about end 1: end moments + shear couple + M0 = 0
        @test q[6] + q[12] + q[8] * L + M0 ≈ 0 atol = 1e-9

        # torsional concentrated moment splits linearly
        qt = fef(N.PointMoment(el, 0.25, [500.0, 0.0, 0.0]))
        @test qt[4] ≈ -500.0 * 0.75 rtol = 1e-12
        @test qt[10] ≈ -500.0 * 0.25 rtol = 1e-12
    end

    @testset "PointMoment in a solve (clamped beam, midspan moment)" begin
        # displacement at the moment location: θ known; test equilibrium and
        # symmetry instead — the two end moment reactions equal M0/4
        nn1 = N.Node([0.0, 0.0, 0.0], :fixed)
        nn2 = N.Node([L, 0.0, 0.0], :fixed)
        e = N.FrameElement(nn1, nn2, sec; rollangle=0.0)
        M0 = 1000.0
        model = N.Model([nn1, nn2], N.AbstractElement{Float64}[e],
            N.AbstractLoad{Float64}[N.PointMoment(e, 0.5, [0.0, 0.0, M0])])
        N.solve!(model)
        R1 = N.reaction(model.results, nn1)
        R2 = N.reaction(model.results, nn2)
        # global moment equilibrium about node 1: M1 + M2 + L·R2y + M0 = 0
        @test R1[6] + R2[6] + L * R2[2] + M0 ≈ 0 atol = 1e-9 * M0
        # antisymmetric shear couple
        @test R1[2] ≈ -R2[2] rtol = 1e-9
    end

    @testset "load case tags" begin
        dead = N.LineLoad(el, [0.0, -2.0, 0.0]; case=:dead)
        live = N.PointLoad(el, 0.5, [0.0, -50.0, 0.0]; case=:live)
        @test dead.case == :dead && live.case == :live
        @test N.NodeForce(n1, [1.0, 0.0, 0.0]).case == :LC1   # default
    end

    @testset "TrapezoidLoad + PointMoment on a VariableElement" begin
        v1 = N.Node([0.0, 0.0, 0.0], :fixed)
        v2 = N.Node([L, 0.0, 0.0], :free)
        vel = N.VariableElement(v1, v2, N.AbstractSection{Float64}[sec, sec], [0.5]; rollangle=0.0)
        # equivalent prismatic model
        p1 = N.Node([0.0, 0.0, 0.0], :fixed)
        p2 = N.Node([L, 0.0, 0.0], :free)
        pel = N.FrameElement(p1, p2, sec; rollangle=0.0)

        for (mkv, mkp) in [
            (N.TrapezoidLoad(vel, 0.2, 0.8, 1.0, 3.0, [0.0, -1.0, 0.0]),
                N.TrapezoidLoad(pel, 0.2, 0.8, 1.0, 3.0, [0.0, -1.0, 0.0])),
            (N.PointMoment(vel, 0.7, [0.0, 0.0, 800.0]),
                N.PointMoment(pel, 0.7, [0.0, 0.0, 800.0])),
        ]
            mv = N.Model([v1, v2], N.AbstractElement{Float64}[vel], N.AbstractLoad{Float64}[mkv])
            mp = N.Model([p1, p2], N.AbstractElement{Float64}[pel], N.AbstractLoad{Float64}[mkp])
            N.solve!(mv; reprocess=true)
            N.solve!(mp; reprocess=true)
            @test N.displacement(mv.results, v2) ≈ N.displacement(mp.results, p2) rtol = 1e-9
            @test N.reaction(mv.results, v1) ≈ N.reaction(mp.results, p1) rtol = 1e-9
        end
    end
end
