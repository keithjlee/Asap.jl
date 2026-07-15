# VariableElement: the super-element must be indistinguishable from (a) a
# single prismatic element when its segments are uniform, and (b) an
# explicit chain of primitive elements with real interior nodes — the
# equivalent system its internal DOFs represent.

@testset "VariableElement" begin
    N = AsapNext
    mat = N.Material(200.0, 77.0, 1.0, 0.3)
    sec = N.Section(mat, 1e4, 8e7, 3e7, 5e6)
    secbig = N.Section(mat, 2e4, 4e8, 1e8, 2e7)

    @testset "construction & geometry" begin
        n1 = N.Node([0.0, 0.0, 0.0], :fixed)
        n2 = N.Node([6000.0, 0.0, 0.0], :free)
        el = N.VariableElement(n1, n2, N.AbstractSection{Float64}[secbig, sec], [0.4], :haunch)

        @test N.n_segments(el) == 2
        @test N.n_internal_dofs(el) == 6
        @test N.ndofs(el) == 18
        @test N.segment_fractions(el) == [0.0, 0.4, 1.0]
        @test N.segment_slots(el, 1) == [collect(1:6); collect(13:18)]
        @test N.segment_slots(el, 2) == [collect(13:18); collect(7:12)]

        @test N.locate_segment(el, 0.2) == (1, 0.5)
        s, τ = N.locate_segment(el, 0.7)
        @test s == 2 && τ ≈ 0.5

        # equal-segment convenience
        el3 = N.VariableElement(n1, n2, N.AbstractSection{Float64}[sec, sec, sec])
        @test N.segment_fractions(el3) ≈ [0.0, 1 / 3, 2 / 3, 1.0]

        @test_throws AssertionError N.VariableElement(n1, n2,
            N.AbstractSection{Float64}[sec, sec], [1.2])

        plain = sprint(show, MIME"text/plain"(), el)
        @test occursin("segments", plain) && occursin("internal DOFs", plain)
    end

    # a two-segment uniform VariableElement must reproduce the single
    # prismatic element exactly (interior nodal values are exact for E-B)
    @testset "uniform segments ≡ prismatic element" begin
        function tip_result(variable::Bool, loadmaker)
            n1 = N.Node([0.0, 0.0, 0.0], :fixed)
            n2 = N.Node([6000.0, 0.0, 0.0], :free)
            el = variable ?
                 N.VariableElement(n1, n2, N.AbstractSection{Float64}[sec, sec], [0.5]) :
                 N.FrameElement(n1, n2, sec)
            model = N.Model([n1, n2], N.AbstractElement{Float64}[el], loadmaker(el))
            N.solve!(model)
            return N.displacement(model.results, n2), N.reaction(model.results, n1), el, model
        end

        cases = [
            el -> N.AbstractLoad{Float64}[N.NodeForce(el.nodeEnd, [0.0, -50.0, 30.0])],
            el -> N.AbstractLoad{Float64}[N.LineLoad(el, [0.0, -2.0, 1.0])],
            el -> N.AbstractLoad{Float64}[N.PointLoad(el, 0.3, [10.0, -40.0, 0.0])],
            el -> N.AbstractLoad{Float64}[N.PointLoad(el, 0.75, [0.0, 0.0, -25.0]),
                N.LineLoad(el, [0.0, -1.0, 0.0])],
        ]
        for loadmaker in cases
            u_v, R_v, elv, mv = tip_result(true, loadmaker)
            u_p, R_p, elp, mp = tip_result(false, loadmaker)
            @test u_v ≈ u_p rtol = 1e-9
            @test R_v ≈ R_p rtol = 1e-9

            # segment end forces at the outer ends match the prismatic element's
            fv = N.element_forces(mv.results, elv, 1)
            fp = N.element_forces(mp.results, elp)
            # atol covers k·u + q̃ cancellation noise at zero-force free ends
            # (force scale here is ~1e5)
            @test fv[1:6] ≈ fp[1:6] rtol = 1e-8 atol = 1e-6
            fv2 = N.element_forces(mv.results, elv, 2)
            @test fv2[7:12] ≈ fp[7:12] rtol = 1e-8 atol = 1e-6
        end
    end

    # a stepped VariableElement must equal the same structure modeled with
    # explicit primitive elements and a real interior node
    @testset "stepped section ≡ explicit element chain" begin
        L, brk = 6000.0, 0.4

        # explicit chain
        m1 = N.Node([0.0, 0.0, 0.0], :fixed)
        mid = N.Node([brk * L, 0.0, 0.0], :free)
        m2 = N.Node([L, 0.0, 0.0], :free)
        e1 = N.FrameElement(m1, mid, secbig)
        e2 = N.FrameElement(mid, m2, sec)
        chain = N.Model([m1, mid, m2], N.AbstractElement{Float64}[e1, e2],
            N.AbstractLoad{Float64}[
                N.NodeForce(m2, [0.0, -50.0, 20.0]),
                N.LineLoad(e1, [0.0, -2.0, 0.0]),      # load on first span only
                N.PointLoad(e2, 0.5, [0.0, 0.0, -30.0]),
            ])
        N.solve!(chain)

        # super-element with the same loads expressed on the whole member:
        # the LineLoad covers [0, 0.4], the point load sits at 0.4 + 0.6·0.5 = 0.7
        v1 = N.Node([0.0, 0.0, 0.0], :fixed)
        v2 = N.Node([L, 0.0, 0.0], :free)
        vel = N.VariableElement(v1, v2, N.AbstractSection{Float64}[secbig, sec], [brk])
        vmodel = N.Model([v1, v2], N.AbstractElement{Float64}[vel],
            N.AbstractLoad{Float64}[
                N.NodeForce(v2, [0.0, -50.0, 20.0]),
                N.DistributedLoad(vel, [0.0, brk], [2.0, 2.0], [0.0, -1.0, 0.0]),
                N.PointLoad(vel, 0.7, [0.0, 0.0, -30.0]),
            ])
        N.solve!(vmodel)

        # tip displacements match; interior-joint displacement (internal DOFs)
        # matches the explicit interior node
        @test N.displacement(vmodel.results, v2) ≈ N.displacement(chain.results, m2) rtol = 1e-9
        part = vmodel.cache.partition
        @test part.n_global == 12 + 6                    # 2 nodes + 1 interior joint
        ioff = vel.internal_offset
        u_internal = vmodel.results.u[ioff+1:ioff+6]
        @test u_internal ≈ collect(N.displacement(chain.results, mid)) rtol = 1e-9

        # reactions and per-segment forces match the explicit chain
        @test N.reaction(vmodel.results, v1) ≈ N.reaction(chain.results, m1) rtol = 1e-9
        @test N.element_forces(vmodel.results, vel, 1) ≈ N.element_forces(chain.results, e1) rtol = 1e-7 atol = 1e-7
        @test N.element_forces(vmodel.results, vel, 2) ≈ N.element_forces(chain.results, e2) rtol = 1e-7 atol = 1e-7
        @test N.axial_force(vmodel.results, vel) ≈ N.axial_force(chain.results, e2) atol = 1e-8
    end

    @testset "self-weight varies with segment area" begin
        L = 6000.0
        n1 = N.Node([0.0, 0.0, 0.0], :fixed)
        n2 = N.Node([L, 0.0, 0.0], :fixed)
        heavy = N.Section(N.Material(200.0, 77.0, 8e-9, 0.3), 2e4, 4e8, 1e8, 2e7)
        light = N.Section(N.Material(200.0, 77.0, 8e-9, 0.3), 1e4, 8e7, 3e7, 5e6)
        vel = N.VariableElement(n1, n2, N.AbstractSection{Float64}[heavy, light], [0.5])
        model = N.Model([n1, n2], N.AbstractElement{Float64}[vel],
            N.AbstractLoad{Float64}[N.SelfWeight(vel; g=[0.0, 0.0, -9810.0])])
        N.solve!(model)

        # total vertical reaction = total weight = Σ ρA·g·Lseg
        W = (N.ρA(heavy) + N.ρA(light)) * 9810.0 * L / 2
        Rz = N.reaction(model.results, n1)[3] + N.reaction(model.results, n2)[3]
        @test Rz ≈ W rtol = 1e-9
        # the heavy half draws more reaction at its support
        @test N.reaction(model.results, n1)[3] > N.reaction(model.results, n2)[3]
    end

    @testset "release on outer end of a super-element" begin
        L = 6000.0
        n1 = N.Node([0.0, 0.0, 0.0], :fixed)
        n2 = N.Node([L, 0.0, 0.0], :pinned)
        vel = N.VariableElement(n1, n2, N.AbstractSection{Float64}[secbig, sec], [0.5];
            release=:fixedfree, rollangle=0.0)
        model = N.Model([n1, n2], N.AbstractElement{Float64}[vel],
            N.AbstractLoad{Float64}[N.LineLoad(vel, [0.0, -2.0, 0.0])])
        N.solve!(model)
        # far end released + pinned support: no moment reaction anywhere at n2
        @test all(abs.(N.reaction(model.results, n2)[4:6]) .< 1e-8)
        # fixed end carries a moment
        @test abs(N.reaction(model.results, n1)[6]) > 1e3
    end
end
