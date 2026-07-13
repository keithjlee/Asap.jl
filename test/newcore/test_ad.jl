# AD validation of the pure functional path: parity with the in-place
# assembler, and Zygote gradients against central finite differences for
# geometry (node position) and section (area) design variables — the
# AsapOptim/DemandTransport use cases, now served by the core.

using Zygote
using FiniteDifferences

const FDM = central_fdm(5, 1)

# AD-test model builders (self-contained, no fixture coupling) --------------

function nc_portal_frame_ad()
    N = AsapNext
    mat = N.Material(200.0, 77.0, 1.0, 0.3)
    sec = N.Section(mat, 1e4, 8e7, 3e7, 5e6)
    n1 = N.Node([0.0, 0.0, 0.0], :fixed)
    n2 = N.Node([0.0, 0.0, 3000.0], :free)
    n3 = N.Node([5000.0, 0.0, 3000.0], :free)
    n4 = N.Node([5000.0, 0.0, 0.0], :fixed)
    beam = N.FrameElement(n2, n3, sec, :beam)
    els = N.AbstractElement{Float64}[
        N.FrameElement(n1, n2, sec, :col), beam, N.FrameElement(n4, n3, sec, :col)]
    return N.Model([n1, n2, n3, n4], els, N.AbstractLoad{Float64}[
        N.LineLoad(beam, [0.0, 0.0, -2.0]),
        N.NodeForce(n2, [5.0, 0.0, 0.0]),
    ])
end

function nc_mixed_ad()
    N = AsapNext
    mat = N.Material(200.0, 77.0, 1.0, 0.3)
    fsec = N.Section(mat, 1e4, 8e7, 3e7, 5e6)
    tsec = N.Section(mat, 1e3)
    n1 = N.Node([0.0, 0.0, 0.0], :fixed)
    n2 = N.Node([3000.0, 0.0, 0.0], :free)
    n3 = N.Node([3000.0, 0.0, 3000.0], :fixed)
    return N.Model([n1, n2, n3],
        N.AbstractElement{Float64}[N.FrameElement(n1, n2, fsec), N.TrussElement(n2, n3, tsec)],
        N.AbstractLoad{Float64}[N.NodeForce(n2, [0.0, 0.0, -100.0])])
end

function nc_truss_ad()
    N = AsapNext
    mat = N.Material(70.0, 1.0, 1.0, 0.3)
    sec = N.Section(mat, 4e-3)
    rot = [true, true, true]
    n1 = N.Node([0.0, 0.0, 0.0], vcat([false, false, false], rot))
    n2 = N.Node([10.0, 0.0, 0.0], vcat([false, false, false], rot))
    n3 = N.Node([5.0, 8.0, 0.0], vcat([true, true, false], rot))
    n4 = N.Node([15.0, 8.0, 0.0], vcat([true, true, false], rot))
    els = N.AbstractElement{Float64}[
        N.TrussElement(n1, n3, sec), N.TrussElement(n2, n3, sec),
        N.TrussElement(n2, n4, sec), N.TrussElement(n3, n4, sec),
    ]
    return N.Model([n1, n2, n3, n4], els, N.AbstractLoad{Float64}[
        N.NodeForce(n3, [0.0, -400.0, 0.0]), N.NodeForce(n4, [200.0, -400.0, 0.0])])
end


@testset "Pure path & AD" begin
    N = AsapNext

    @testset "parity: pure assemble_K ≡ in-place assemble_K!" begin
        for build in (nc_portal_frame_ad, nc_mixed_ad)
            model = build()
            N.process!(model)
            cache = model.cache
            N.assemble_loads!(cache, model)

            K_inplace = copy(N.assemble_K!(cache, model))
            K_pure = N.assemble_K(cache, N.extract_state(model))
            @test K_pure ≈ K_inplace rtol = 1e-13

            # pure solve ≡ in-place solve
            N.solve!(model)
            u_pure = N.solve(model, N.extract_state(model))
            @test u_pure ≈ model.results.u rtol = 1e-11
        end
    end

    @testset "gradients w.r.t. node positions (Zygote vs FiniteDifferences)" begin
        model = nc_portal_frame_ad()
        N.solve!(model)                                # builds cache + loads
        state0 = N.extract_state(model)
        X0 = vec(state0.X)

        function obj(Xflat)
            return N.compliance(model, N.ModelState{Float64}(reshape(Xflat, 3, :), state0.sections))
        end

        g_zygote = Zygote.gradient(obj, X0)[1]
        g_fd = FiniteDifferences.grad(FDM, obj, X0)[1]
        @test g_zygote ≈ g_fd rtol = 1e-6
        @test norm(g_zygote) > 0
    end

    @testset "gradients w.r.t. section areas (truss sizing)" begin
        model = nc_truss_ad()
        N.solve!(model)
        state0 = N.extract_state(model)
        mat = N.Material(70.0, 1.0, 1.0, 0.3)
        A0 = fill(4e-3, length(model.elements))

        function obj(A)
            sections = map(a -> N.Section(mat, a), A)   # map, not a typed
            # comprehension — typed collection fills by setindex! (Zygote-hostile)
            return N.compliance(model, N.ModelState{Float64}(state0.X, sections))
        end

        g_zygote = Zygote.gradient(obj, A0)[1]
        g_fd = FiniteDifferences.grad(FDM, obj, A0)[1]
        @test g_zygote ≈ g_fd rtol = 1e-6
        # compliance decreases with added material: all sensitivities negative
        @test all(<(0), g_zygote)
    end

    @testset "gradients w.r.t. connection stiffness (semi-rigid design)" begin
        # differentiable joint stiffness — the Forma concrete use case
        mat = N.Material(200.0, 77.0, 1.0, 0.3)
        sec = N.Section(mat, 1e4, 8e7, 3e7, 5e6)
        n1 = N.Node([0.0, 0.0, 0.0], :fixed)
        n2 = N.Node([3000.0, 0.0, 0.0], :free)

        function tip_deflection(kz)
            # Ψ = 0 so the global-Y load bends in the local x–y plane (EIx),
            # governed by the kz springs
            ends = N.EndConditions(N.EndSprings(Inf, Inf, Inf, kz), N.rigid_end())
            el = N.FrameElement(n1, n2, sec, ends, :beam; Ψ=0.0)
            model = N.Model([n1, n2], N.AbstractElement{Float64}[el],
                N.AbstractLoad{Float64}[N.NodeForce(n2, [0.0, -10.0, 0.0])])
            N.solve!(model)
            return N.displacement(model.results, n2)[2]
        end

        # forward physics: stiffer base connection → smaller tip deflection,
        # bracketed by the ideal-release and rigid limits
        δ_soft = abs(tip_deflection(1e5))
        δ_stiff = abs(tip_deflection(1e8))
        δ_rigid = abs(tip_deflection(Inf))
        L, P = 3000.0, 10.0
        @test δ_rigid ≈ P * L^3 / (3 * N.EIx(sec)) rtol = 1e-9
        @test δ_rigid < δ_stiff < δ_soft
        # base-spring analytic solution: δ = PL³/3EI + (PL)·L/k
        @test δ_stiff ≈ P * L^3 / (3 * N.EIx(sec)) + P * L * L / 1e8 rtol = 1e-9

        # spring-as-design-variable flows through fixity_factor smoothly:
        p(kz) = N.fixity_factor(kz, N.EIx(sec), 3000.0)
        # rtol reflects finite-difference truncation at the 1e8 argument scale
        @test Zygote.gradient(p, 1e8)[1] ≈ FiniteDifferences.grad(FDM, p, 1e8)[1] rtol = 1e-4
    end
end

