# Shell S6b: feature parity — plot extraction (ShellGeo, ModelGeo on mixed
# models), principal forces, shell envelopes, element connectivity, display.

using StaticArrays, LinearAlgebra

@testset "Shell feature parity (S6b)" begin
    conc = Material(3.0e7, 2.4, 0.2)
    w = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 4.0, 3.0, ShellSection(conc, 0.2); size = 1.0)
    top = [w.grid[i, end] for i in eachindex(w.s)]
    tip = Node([6.0, 0.0, 3.0], :free)
    beam = FrameElement(top[end], tip, Section(Material(2.0e8, 7.85, 0.3), 1e-2, 1e-4, 1e-4, 1e-5), :beam)
    els = AbstractElement{Float64}[w.elements..., beam]
    loads = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81], case = :D) for el in w.elements]...,
        NodeForce(tip, [0.0, 0, -5.0]; case = :D), [NodeForce(n, [4.0, 0, 0]; case = :W) for n in top]...]
    m = Model(vcat(w.nodes, tip), els, loads)
    solve!(m)

    @testset "ShellGeo and nodal averaging" begin
        g = ShellGeo(m)
        @test size(g.quads) == (4, length(w.elements))
        @test g.elements == collect(1:length(w.elements))
        @test g.N[:, 1] ≈ shell_resultants(m, w.elements[1]).N
        @test all(isapprox.(g.normals[2, :], -1.0))
        @test all(g.thickness .== 0.2)
        avg = nodal_average(g, g.N[2, :])
        @test length(avg) == length(m.nodes) && avg[end] == 0.0       # the beam tip: no shell
        @test_throws DimensionMismatch nodal_average(g, [1.0])
        @test size(g.disp) == (3, length(m.nodes))
    end

    @testset "principal membrane forces" begin
        N1, N2, θ = principal_forces([3.0, -1.0, 2.0])
        R = [cos(θ) sin(θ); -sin(θ) cos(θ)]
        @test R * [3.0 2.0; 2.0 -1.0] * R' ≈ [N1 0; 0 N2] atol = 1e-12
        @test N1 >= N2
    end

    @testset "ModelGeo on a mixed model covers the line elements" begin
        mg = Geo(m)
        @test length(mg.indices) == 1 && mg.indices[1] == [beam.nodeStart.index, beam.nodeEnd.index]
        @test mg.lengths ≈ [2.0]
    end

    @testset "shell envelopes over combinations" begin
        cr = solve_cases!(m)
        combos = [LoadCombination(:A, [:D => 1.2, :W => 1.6]), LoadCombination(:B, [:D => 0.9, :W => -1.0])]
        el = w.elements[2]
        env = shell_envelope(m, el, cr, combos; at = (0.5, -0.5))
        rA = shell_resultants(m, el; at = (0.5, -0.5), results = combine(cr, combos[1]), factors = Dict(:D => 1.2, :W => 1.6))
        rB = shell_resultants(m, el; at = (0.5, -0.5), results = combine(cr, combos[2]), factors = Dict(:D => 0.9, :W => -1.0))
        @test env.Nmax ≈ max.(rA.N, rB.N) && env.Nmin ≈ min.(rA.N, rB.N)
        @test env.Mmax ≈ max.(rA.M, rB.M) && env.Qmin ≈ min.(rA.Q, rB.Q)
        @test_throws ArgumentError shell_envelope(m, el, cr, LoadCombination{Float64}[])
    end

    @testset "element connectivity and display" begin
        C = element_connectivity(m)
        b = length(els)
        @test C[b, findfirst(el -> any(n -> n === top[end], el.nodes), w.elements)] == 1
        @test C[1, 1] == 0 && all(C.nzval .== 1)
        str = sprint(show, MIME"text/plain"(), m)
        @test occursin("shell", str)
    end
    @testset "nodal_average weights: count, area, angle" begin
        # a 2 × 1 quad (area 2) and a right triangle (area 1/2) sharing nodes (2,0) and (2,1)
        P = [0.0 2.0 2.0 0.0 3.0; 0.0 0.0 1.0 1.0 0.0]
        sec = ShellSection(Material(3.0e4, 2.4e-3, 0.2), 0.2)
        mesh = shell_mesh(P, [(1, 2, 3, 4), (2, 5, 3)], sec)
        foreach(n -> fixnode!(n, :fixed), mesh.nodes[[1, 4]])
        m = Model(mesh.nodes, AbstractElement{Float64}[mesh.elements...],
            AbstractLoad{Float64}[NodeForce(mesh.nodes[5], [1.0, 0, 0])])
        solve!(m)
        g = ShellGeo(m)
        v = [el isa ShellElement{Float64,4} ? 1.0 : 3.0 for el in m.elements[g.elements]]
        n3 = findfirst(n -> n.position ≈ [2.0, 1, 0], m.nodes)   # quad corner 90°, triangle corner 45°
        n1 = findfirst(n -> n.position ≈ [0.0, 0, 0], m.nodes)   # quad only
        @test nodal_average(g, v)[n3] ≈ 2.0
        @test nodal_average(g, v; weight = :area)[n3] ≈ (2 * 1 + 0.5 * 3) / 2.5
        @test nodal_average(g, v; weight = :angle)[n3] ≈ (π / 2 * 1 + π / 4 * 3) / (3π / 4)
        @test all(w -> nodal_average(g, v; weight = w)[n1] ≈ 1.0, (:count, :area, :angle))
        @test_throws ArgumentError nodal_average(g, v; weight = :volume)
    end
end