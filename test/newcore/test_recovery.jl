# Phase 3: equilibrium internal-force recovery.
#
# Validation layers:
#   1. closed-form beam diagrams (cantilever, simply supported)
#   2. calculus property tests: dMz/dx = Vy, dVy/dx = wy_local, dMy/dx = −Vz,
#      and exact jumps at point actions
#   3. the pinned AsapToolkit InternalForces oracles on the portal frames,
#      under the documented axis-name map (tk.My→Mz, tk.Mz→−My, tk.P→N)
#   4. displacement recovery against classic deflection formulas
#   5. VariableElement continuity and equivalence with explicit chains

@testset "Force recovery (Phase 3)" begin
    N = AsapNext
    mat = N.Material(200.0, 77.0, 1.0, 0.3)
    sec = N.Section(mat, 1e4, 8e7, 3e7, 5e6)
    L = 4000.0

    function cantilever(loads_of; rollangle=0.0)
        n1 = N.Node([0.0, 0.0, 0.0], :fixed)
        n2 = N.Node([L, 0.0, 0.0], :free)
        el = N.FrameElement(n1, n2, sec; rollangle=rollangle)
        model = N.Model([n1, n2], N.AbstractElement{Float64}[el], loads_of(el, n2))
        N.solve!(model)
        return model, el
    end

    @testset "cantilever diagrams (analytic)" begin
        P = 100.0
        model, el = cantilever((el, tip) -> N.AbstractLoad{Float64}[
            N.NodeForce(tip, [0.0, -P, 0.0])])
        st = N.internal_forces(model, el)

        for t in (0.0, 0.25, 0.5, 0.9)
            @test N.shear_y(st, t) ≈ P rtol = 1e-9           # constant shear
            @test N.moment_z(st, t) ≈ -P * L * (1 - t) rtol = 1e-9  # hogging
        end
        @test abs(N.moment_z(st, 1.0)) < 1e-8 * P * L         # free end

        # UDL cantilever: Mz(x) = −w(L−x)²/2
        w = 2.0
        model2, el2 = cantilever((el, tip) -> N.AbstractLoad{Float64}[
            N.LineLoad(el, [0.0, -w, 0.0])])
        st2 = N.internal_forces(model2, el2)
        for t in (0.0, 0.3, 0.7, 1.0)
            x = t * L
            @test N.moment_z(st2, t) ≈ -w * (L - x)^2 / 2 rtol = 1e-9 atol = 1e-8
            @test N.shear_y(st2, t) ≈ w * (L - x) rtol = 1e-9 atol = 1e-8
        end
    end

    @testset "simply supported UDL: midspan Mz = wL²/8" begin
        w = 3.0
        n1 = N.Node([0.0, 0.0, 0.0], :pinned)
        n2 = N.Node([L, 0.0, 0.0], [true, false, false, true, true, true])
        el = N.FrameElement(n1, n2, sec; rollangle=0.0)
        model = N.Model([n1, n2], N.AbstractElement{Float64}[el],
            N.AbstractLoad{Float64}[N.LineLoad(el, [0.0, -w, 0.0])])
        N.planarize!(model)
        N.solve!(model)
        st = N.internal_forces(model, el)
        @test N.moment_z(st, 0.5) ≈ w * L^2 / 8 rtol = 1e-9            # sagging +
        @test abs(N.moment_z(st, 0.0)) < 1e-8 * w * L^2
        @test abs(N.moment_z(st, 1.0)) < 1e-8 * w * L^2
        @test N.shear_y(st, 0.0) ≈ w * L / 2 rtol = 1e-9
    end

    @testset "calculus properties + jumps" begin
        # rich loading: trapezoid + point force + point moment, skewed element
        n1 = N.Node([0.0, 0.0, 0.0], :fixed)
        n2 = N.Node([3000.0, 1000.0, 2000.0], :fixed)
        el = N.FrameElement(n1, n2, sec)
        Lel = Base.length(el)
        loads = N.AbstractLoad{Float64}[
            N.TrapezoidLoad(el, 0.1, 0.8, 1.0, 4.0, normalize([0.2, -0.9, 0.4])),
            N.PointLoad(el, 0.45, [50.0, -80.0, 30.0]),
            N.PointMoment(el, 0.65, [2e4, -1e4, 3e4]),
        ]
        model = N.Model([n1, n2], N.AbstractElement{Float64}[el], loads)
        N.solve!(model)
        st = N.internal_forces(model, el)

        # derivative identities via central differences at smooth points
        h = 1e-6
        for t in (0.2, 0.3, 0.55, 0.9)
            x = t * Lel
            dMz = (N.moment_z(st, t + h) - N.moment_z(st, t - h)) / (2h * Lel)
            @test dMz ≈ N.shear_y(st, t) rtol = 1e-5
            dMy = (N.moment_y(st, t + h) - N.moment_y(st, t - h)) / (2h * Lel)
            @test dMy ≈ -N.shear_z(st, t) rtol = 1e-5
            dVy = (N.shear_y(st, t + h) - N.shear_y(st, t - h)) / (2h * Lel)
            @test dVy ≈ N._trace_w(st.trace, 2, x) rtol = 1e-4 atol = 1e-10
        end

        # exact jumps at the point force (local components)
        Λ = N.local_frame(el)
        Ploc = Λ * SVector(50.0, -80.0, 30.0)
        ε = 1e-10
        @test N.shear_y(st, 0.45 + ε) - N.shear_y(st, 0.45 - ε) ≈ Ploc[2] rtol = 1e-6
        @test N.shear_z(st, 0.45 + ε) - N.shear_z(st, 0.45 - ε) ≈ Ploc[3] rtol = 1e-6
        @test N.axial_force(st, 0.45 + ε) - N.axial_force(st, 0.45 - ε) ≈ -Ploc[1] rtol = 1e-6

        # exact jumps at the point moment
        Mloc = Λ * SVector(2e4, -1e4, 3e4)
        @test N.moment_z(st, 0.65 + ε) - N.moment_z(st, 0.65 - ε) ≈ -Mloc[3] rtol = 1e-6
        @test N.torsion(st, 0.65 + ε) - N.torsion(st, 0.65 - ε) ≈ -Mloc[1] rtol = 1e-6

        # end-value ties: recovered forces at t=1 equal the stored end block
        f = N.element_forces(model.results, el)
        @test N.axial_force(st, 1.0) ≈ f[7] rtol = 1e-8
        @test N.torsion(st, 1.0) ≈ f[10] rtol = 1e-8
    end

    @testset "vs pinned AsapToolkit oracles (portal frames)" begin
        # documented axis-name map: tk.P → N, tk.Vy → Vy, tk.My → Mz,
        # tk.Vz → Vz, tk.Mz → −My
        for r in CHAR_RELEASES
            model = nc_portal_frame(r)
            N.solve!(model)
            for el in model.elements
                tag = "portal_$r/$(el.id)$(el.index)"
                haskey(TOOLKIT_FIXTURES, "$tag/x") || continue   # joist beam gap
                xs = TOOLKIT_FIXTURES["$tag/x"]
                st = N.internal_forces(model, el)
                Lel = st.L
                # skip stations at point-action jumps (left/right ambiguity)
                jumps = st.trace.pstation
                keep = [i for (i, x) in enumerate(xs) if all(abs(x - a) > 1e-6 * Lel for a in jumps)]
                ts = xs[keep] ./ Lel

                @test [N.axial_force(st, t) for t in ts] ≈ TOOLKIT_FIXTURES["$tag/P"][keep] rtol = 1e-6 atol = 1e-6
                @test [N.shear_y(st, t) for t in ts] ≈ TOOLKIT_FIXTURES["$tag/Vy"][keep] rtol = 1e-6 atol = 1e-6
                @test [N.moment_z(st, t) for t in ts] ≈ TOOLKIT_FIXTURES["$tag/My"][keep] rtol = 1e-6 atol = 1e-5
                @test [N.shear_z(st, t) for t in ts] ≈ TOOLKIT_FIXTURES["$tag/Vz"][keep] rtol = 1e-6 atol = 1e-6
                @test [-N.moment_y(st, t) for t in ts] ≈ TOOLKIT_FIXTURES["$tag/Mz"][keep] rtol = 1e-6 atol = 1e-5
            end
        end
    end

    @testset "displacement recovery (classic deflections)" begin
        # cantilever tip load: v(L) = −PL³/3EI, midspan v = −5PL³/48EI
        P = 100.0
        model, el = cantilever((el, tip) -> N.AbstractLoad{Float64}[
            N.NodeForce(tip, [0.0, -P, 0.0])])
        st = N.internal_forces(model, el)
        EIx = N.EIx(sec)
        @test N.local_displacements(st, 1.0)[2] ≈ -P * L^3 / (3EIx) rtol = 1e-8
        @test N.local_displacements(st, 0.5)[2] ≈ -5P * L^3 / (48EIx) rtol = 1e-8

        # recovery at t=1 matches the nodal solution
        tip_u = N.displacement(model.results, el.nodeEnd)
        @test N.local_displacements(st, 1.0)[2] ≈ tip_u[2] rtol = 1e-9

        # simply supported UDL: midspan v = −5wL⁴/384EI
        w = 3.0
        n1 = N.Node([0.0, 0.0, 0.0], :pinned)
        n2 = N.Node([L, 0.0, 0.0], [true, false, false, true, true, true])
        e2 = N.FrameElement(n1, n2, sec; rollangle=0.0)
        m2 = N.Model([n1, n2], N.AbstractElement{Float64}[e2],
            N.AbstractLoad{Float64}[N.LineLoad(e2, [0.0, -w, 0.0])])
        N.planarize!(m2)
        N.solve!(m2)
        st2 = N.internal_forces(m2, e2)
        @test N.local_displacements(st2, 0.5)[2] ≈ -5w * L^4 / (384EIx) rtol = 1e-8
    end

    @testset "displacement recovery at RELEASED/semi-rigid start ends" begin
        #=
        The member-side start rotation differs from the NODE rotation by the
        hinge/spring jump — the recovery derives it from far-end
        compatibility. Regression for the :joist case: node rotations fixed,
        member internally simply supported.
        =#
        w = 3.0
        EIx = N.EIx(sec)

        # :joist between rotation-restrained nodes = simply supported inside
        j1 = N.Node([0.0, 0.0, 0.0], :fixed)
        j2 = N.Node([L, 0.0, 0.0], :fixed)
        ej = N.FrameElement(j1, j2, sec, N.EndConditions(:joist), :beam; rollangle=0.0)
        mj = N.Model([j1, j2], N.AbstractElement{Float64}[ej],
            N.AbstractLoad{Float64}[N.LineLoad(ej, [0.0, -w, 0.0])])
        N.solve!(mj)
        stj = N.internal_forces(mj, ej)
        @test N.local_displacements(stj, 0.5)[2] ≈ -5w * L^4 / (384EIx) rtol = 1e-8
        # far-end compatibility: recovered curve lands on the end node
        @test abs(N.local_displacements(stj, 1.0)[2]) < 1e-10

        # :freefixed (hinge at the START): propped cantilever,
        # midspan v = −wL⁴/192EI
        p1 = N.Node([0.0, 0.0, 0.0], :fixed)
        p2 = N.Node([L, 0.0, 0.0], :fixed)
        ep = N.FrameElement(p1, p2, sec, N.EndConditions(:freefixed), :beam; rollangle=0.0)
        mp = N.Model([p1, p2], N.AbstractElement{Float64}[ep],
            N.AbstractLoad{Float64}[N.LineLoad(ep, [0.0, -w, 0.0])])
        N.solve!(mp)
        stp = N.internal_forces(mp, ep)
        @test N.local_displacements(stp, 0.5)[2] ≈ -w * L^4 / (192EIx) rtol = 1e-8

        # semi-rigid springs sweep: deflection decreases monotonically with k
        # and spans the released → rigid limits
        function δmid(k)
            s1 = N.Node([0.0, 0.0, 0.0], :fixed)
            s2 = N.Node([L, 0.0, 0.0], :fixed)
            ends = N.EndConditions(N.EndSprings(Inf, Inf, k, k), N.EndSprings(Inf, Inf, k, k))
            es = N.FrameElement(s1, s2, sec, ends, :beam; rollangle=0.0)
            ms = N.Model([s1, s2], N.AbstractElement{Float64}[es],
                N.AbstractLoad{Float64}[N.LineLoad(es, [0.0, -w, 0.0])])
            N.solve!(ms)
            abs(N.local_displacements(N.internal_forces(ms, es), 0.5)[2])
        end
        δs = δmid.([1e2, 1e4, 1e6, 1e8])
        @test issorted(δs; rev=true)
        @test δs[1] < 5w * L^4 / (384EIx) * (1 + 1e-6)
        @test δs[end] > w * L^4 / (384EIx) * (1 - 1e-6)
    end

    @testset "VariableElement: unified queries + interior continuity" begin
        n1 = N.Node([0.0, 0.0, 0.0], :fixed)
        n2 = N.Node([L, 0.0, 0.0], :free)
        big = N.Section(mat, 2e4, 4e8, 1e8, 2e7)
        vel = N.VariableElement(n1, n2, N.AbstractSection{Float64}[big, sec], [0.4]; rollangle=0.0)
        model = N.Model([n1, n2], N.AbstractElement{Float64}[vel],
            N.AbstractLoad{Float64}[
                N.LineLoad(vel, [0.0, -2.0, 0.0]),
                N.PointLoad(vel, 0.7, [0.0, -60.0, 0.0]),
            ])
        N.solve!(model)

        # moment and shear are continuous across the interior joint (no load there)
        @test N.moment_z(model, vel, 0.4 - 1e-9) ≈ N.moment_z(model, vel, 0.4 + 1e-9) rtol = 1e-6
        @test N.shear_y(model, vel, 0.4 - 1e-9) ≈ N.shear_y(model, vel, 0.4 + 1e-9) rtol = 1e-6

        # unified fraction queries match the analytic cantilever solution:
        # Mz(x) = −[w(L−x)²/2 + P·(xP−x)⁺]
        w, P, xP = 2.0, 60.0, 0.7L
        for t in (0.1, 0.4, 0.55, 0.85)
            x = t * L
            expected = -(w * (L - x)^2 / 2 + (x < xP ? P * (xP - x) : 0.0))
            @test N.moment_z(model, vel, t) ≈ expected rtol = 1e-8
        end

        # dense sampling spans segments with increasing x
        f = N.InternalForces(model, vel; resolution=24)
        @test issorted(f.x)
        @test first(f.x) == 0.0 && last(f.x) ≈ L
    end

    @testset "scalar evaluators allocation-free" begin
        model, el = cantilever((el, tip) -> N.AbstractLoad{Float64}[
            N.LineLoad(el, [0.0, -2.0, 0.0])])
        st = N.internal_forces(model, el)
        alloc(f, s, t) = (f(s, t); @allocated f(s, t))
        @test alloc(N.moment_z, st, 0.37) == 0
        @test alloc(N.shear_y, st, 0.37) == 0
        @test alloc(N.axial_force, st, 0.37) == 0
    end
end
