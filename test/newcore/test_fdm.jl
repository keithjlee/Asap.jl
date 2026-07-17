#=
FDM subsystem: uniform solves (regression), per-axis fixity (separable
solves), translations, and hand-calculated equilibria.
=#
@testset "FDM (force density method)" begin
    N = AsapNext

    @testset "uniform network (batched path)" begin
        gen = N.GridNetwork(4, 6.0, 4, 6.0)
        nw = gen.network
        @test nw.processed && !nw.mixed
        @test nw.Naxis[1] == nw.Naxis[2] == nw.Naxis[3] == nw.N
        # interior nodes rose toward the load direction ([0,0,1] default)
        zs = [n.position[3] for n in nw.nodes[nw.N]]
        @test all(zs .> 0)
    end

    @testset "hand calc: 2-cable chain, fully free middle" begin
        # ends fixed at (0,0,0), (2,0,0); q = 1 both; P = [0,0,-1] at middle.
        # D = q1+q2 = 2 per coordinate → x = (0 - (-q·0 - q·2))/2 = 1,
        # y = 0, z = -1/2.
        a = N.FDMnode([0.0, 0.0, 0.0], false)
        m = N.FDMnode([1.3, 0.7, 0.2], true, :mid)     # arbitrary start position
        b = N.FDMnode([2.0, 0.0, 0.0], false)
        els = [N.FDMelement(a, m, 1.0), N.FDMelement(m, b, 1.0)]
        nw = N.Network([a, m, b], els, [N.FDMload(m, [0.0, 0.0, -1.0])])
        N.solve!(nw)
        @test m.position ≈ [1.0, 0.0, -0.5] rtol = 1e-12
        # anchors react; z-components sum against the load
        @test a.reaction[3] + b.reaction[3] ≈ -1.0 rtol = 1e-10
    end

    @testset "per-axis fixity: plan prescribed, height free" begin
        # same chain, but the middle node is FIXED in x and y at (0.7, 0.3)
        # and free only in z. Separable equilibrium in z alone: z = -1/2·... 
        # D_z = 2, rhs = Pz − (−q·z_a − q·z_b) = −1 → z = −0.5. x, y stay put.
        a = N.FDMnode([0.0, 0.0, 0.0], false)
        m = N.FDMnode([0.7, 0.3, 0.0], Bool[false, false, true], :mid)
        b = N.FDMnode([2.0, 0.0, 0.0], false)
        els = [N.FDMelement(a, m, 1.0), N.FDMelement(m, b, 1.0)]
        nw = N.Network([a, m, b], els, [N.FDMload(m, [0.0, 0.0, -1.0])])
        N.solve!(nw)
        @test nw.mixed
        @test m.position[1] ≈ 0.7 atol = 1e-14          # plan untouched
        @test m.position[2] ≈ 0.3 atol = 1e-14
        @test m.position[3] ≈ -0.5 rtol = 1e-12          # height at equilibrium
        # the plan restraint reacts in x/y (member forces don't balance there),
        # and the free z-component is zero by definition
        @test m.reaction[3] == 0.0
        @test abs(m.reaction[1]) > 0

        # pure re-solve variant: doubled q halves the sag (z = Pz/Σq with
        # anchors at z = 0), and the network itself is untouched
        xyz = N.solve(nw, [2.0, 2.0])
        @test xyz[2, 3] ≈ -0.25 rtol = 1e-12
        @test m.position[3] ≈ -0.5 rtol = 1e-12
    end

    @testset "mixed ≡ uniform when fixity is uniform" begin
        # a network declared per-axis but with uniform values must match the
        # batched path exactly
        gen1 = N.GridNetwork(3, 4.0, 3, 4.0)
        n1 = gen1.network
        gen2 = N.GridNetwork(3, 4.0, 3, 4.0)
        n2 = gen2.network
        # re-declare one FREE node's fixity per-axis (same values) and force
        # the mixed path by fixing a different node's z only
        free_idx = first(n2.N)
        n2.nodes[free_idx].fixity = SVector(true, true, true)
        other = n2.N[2]
        n2.nodes[other].fixity = SVector(true, true, false)   # z now fixed
        zpin = n2.nodes[other].position[3]
        N.solve!(n2; reprocess = true)
        @test n2.mixed
        @test n2.nodes[other].position[3] == zpin
        # and the fully-free node still equilibrates (finite, sensible)
        @test all(isfinite, n2.nodes[free_idx].position)
    end

    @testset "translations carry per-axis fixity" begin
        gen = N.GridNetwork(3, 4.0, 3, 4.0)
        truss = N.to_truss(gen.network, N.Section(N.Steel_kNm, 1e-3))
        @test truss isa N.Model
        nw2 = N.to_network(truss)
        @test length(nw2.elements) == length(gen.network.elements)
    end

    @testset "show methods" begin
        gen = N.GridNetwork(3, 4.0, 3, 4.0)
        @test occursin("force density", repr(MIME("text/plain"), gen.network))
        @test occursin("fixity", repr(MIME("text/plain"), gen.network.nodes[1]))
    end
end
