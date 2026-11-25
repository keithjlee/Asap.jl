# S7e/S7f: general planar meshes — the array constructor (any mesher), polygon meshing
# with holes (DelaunayTriangulation.jl extension), boundary-edge node sets, consistent
# triangle edge loads, segment queries and reactions per bearing on a deep beam.

using LinearAlgebra, StaticArrays
import DelaunayTriangulation

@testset "General shell meshes (S7e, S7f)" begin
    sec = ShellSection(Material(3.0e4, 2.4e-3, 0.2), 0.3)

    @testset "shell_mesh from arrays: mixed cells, orientation, boundary loops" begin
        # 3 × 2 grid of points, two quads on the left, two triangles on the right (one CW)
        P = [0.0 1.0 2.0 0.0 1.0 2.0; 0.0 0.0 0.0 1.0 1.0 1.0]
        cells = [(1, 2, 5, 4), (2, 3, 6), (2, 6, 5)]          # (2, 6, 5) is clockwise
        m = shell_mesh(P, cells, sec; origin = [10.0, 0, 0], e1 = [0.0, 1, 0], e2 = [0.0, 0, 1])
        @test length(m.nodes) == 6
        @test count(el -> el isa ShellElement{Float64,4}, m.elements) == 1
        @test count(el -> el isa ShellElement{Float64,3}, m.elements) == 2
        @test all(el -> local_frame(el)[3, :] ≈ [1.0, 0, 0], m.elements)   # e3 = e1 × e2
        @test m.nodes[3].position ≈ [10.0, 2.0, 0.0]
        loop = boundary_nodes(m)                          # outline, counter-clockwise
        @test length(loop) == 6
        @test Set(loop) == Set(m.nodes)
        @test [n.position for n in nodes_on_segment(m, (0.0, 0.0), (2.0, 0.0))] ≈
              [[10.0, 0, 0], [10.0, 1, 0], [10.0, 2, 0]]
    end

    @testset "mesh_polygon: holes, quality, tracked edges" begin
        outline = [(0.0, 0.0), (6.0, 0.0), (6.0, 3.0), (0.0, 3.0)]
        circle = circle_polygon((1.6, 1.5), 0.6; n = 24)
        L = [(3.5, 0.8), (5.0, 0.8), (5.0, 1.4), (4.1, 1.4), (4.1, 2.2), (3.5, 2.2)]   # L-shaped
        m = mesh_polygon(outline, sec; holes = [circle, L], h = 0.25, min_angle = 25)
        @test all(el -> el isa ShellElement{Float64,3}, m.elements)
        # area of the mesh = outline − holes (circle as its 24-gon)
        Acirc = 24 / 2 * 0.6^2 * sin(2π / 24)
        AL = 1.5 * 0.6 + 0.6 * 0.8
        @test sum(Asap.area, m.elements) ≈ 18.0 - Acirc - AL rtol = 1e-10
        # quality: every angle ≥ 25°
        minang = minimum(m.elements) do el
            x = map(n -> n.position, el.nodes)
            minimum(k -> acosd(clamp(dot(normalize(x[mod1(k + 1, 3)] - x[k]), normalize(x[mod1(k + 2, 3)] - x[k])), -1, 1)), 1:3)
        end
        @test minang >= 25 - 1e-6
        # boundary: 1 outline + 2 holes; outline edge 2 (x = 6) runs (6,0) → (6,3)
        @test length(m.boundary) == 3
        right = boundary_nodes(m, 1, 2)
        @test right[1].position ≈ [6.0, 0, 0] && right[end].position ≈ [6.0, 3, 0]
        @test all(n -> n.position[1] ≈ 6.0, right)
        @test length(right) >= 3                                       # refined
        @test Set(right) == Set(nodes_on_segment(m, (6.0, 0.0), (6.0, 3.0)))
        @test length(boundary_nodes(m, 2)) >= 24                       # circle: all its vertices
        # a hole along the outline edge (a notch) works; touching at a vertex is refused clearly
        notch = mesh_polygon(outline, sec; holes = [[(2.0, 0.0), (3.0, 0.0), (3.0, 1.0), (2.0, 1.0)]], h = 0.5)
        @test sum(Asap.area, notch.elements) ≈ 17.0 rtol = 1e-10
        @test_throws ArgumentError mesh_polygon(outline, sec; holes = [[(0.0, 0.0), (1.0, 0.5), (0.5, 1.0)]], h = 0.5)
        # an extra point becomes a node
        mp = mesh_polygon(outline, sec; h = 0.5, points = [(3.0, 1.5)])
        @test any(n -> n.position ≈ [3.0, 1.5, 0.0], mp.nodes)
    end

    @testset "triangle edge_loads: uniform stress reproduced exactly" begin
        # panel in tension σxx = p/t through edge loads on x = 4; held on x = 0
        outline = [(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (0.0, 2.0)]
        t = 0.3
        m = mesh_polygon(outline, sec; h = 0.6)
        for n in m.nodes
            x, y = n.position[1], n.position[2]
            # held edge: ux (uy at one node) and the drilling rotation, which a uniform
            # state keeps at zero (as for any drilling membrane, a held edge fixes ψ)
            fixnode!(n, x ≈ 0 ? Bool[0, y ≈ 0 ? 0 : 1, 0, 0, 0, 0] : Bool[1, 1, 0, 0, 0, 1])
        end
        p = 5.0                                          # [force/length] along the edge
        loads = AbstractLoad{Float64}[edge_loads(m, (4.0, 0.0), (4.0, 2.0), [p, 0.0, 0.0])...]
        model = Model(m.nodes, AbstractElement{Float64}[m.elements...], loads)
        solve!(model)
        for el in m.elements
            r = shell_resultants(model, el)
            Q = r.frame[1:2, 1:2]
            Ng = Q' * [r.N[1] r.N[3]; r.N[3] r.N[2]] * Q
            @test Ng ≈ [p 0; 0 0] atol = 1e-8
        end
        # reaction total on the held edge balances the load
        F, _ = total_reaction(model.results, nodes_on_segment(m, (0.0, 0.0), (0.0, 2.0)))
        @test F ≈ [-2p, 0, 0] atol = 1e-8
    end

    @testset "deep beam on two bearings: statics of the reactions" begin
        # 6 × 3 panel, bearings 0–0.5 and 5.5–6 on the bottom, line load on the top
        # 2–4, self-weight; vertices at every bearing/load end point
        outline = [(0.0, 0.0), (0.5, 0.0), (5.5, 0.0), (6.0, 0.0), (6.0, 3.0),
                   (4.0, 3.0), (2.0, 3.0), (0.0, 3.0)]
        m = mesh_polygon(outline, sec; h = 0.3, e1 = [1.0, 0, 0], e2 = [0.0, 0, 1])
        # wall in the XZ plane: out-of-plane (Y) translation and in-plane-axis rotations held
        for n in m.nodes
            fixnode!(n, Bool[1, 0, 1, 0, 1, 0])
        end
        brgL = nodes_on_segment(m, (0.0, 0.0), (0.5, 0.0))
        brgR = nodes_on_segment(m, (5.5, 0.0), (6.0, 0.0))
        for n in brgL
            fixnode!(n, Bool[0, 0, 0, 0, 0, 0])            # pinned bearing (held line: ψ too)
        end
        for n in brgR
            fixnode!(n, Bool[1, 0, 0, 0, 0, 0])            # roller in x
        end
        w = 40.0
        loads = AbstractLoad{Float64}[edge_loads(m, (2.0, 3.0), (4.0, 3.0), [0.0, 0.0, -w])...,
            [SelfWeight(el; g = [0.0, 0, -9.81]) for el in m.elements]...]
        model = Model(m.nodes, AbstractElement{Float64}[m.elements...], loads)
        solve!(model)
        W = areal_mass(sec) * 9.81 * 18.0
        FL, _ = total_reaction(model.results, brgL)
        FR, _ = total_reaction(model.results, brgR)
        @test FL[3] + FR[3] ≈ 2w + W rtol = 1e-9
        @test abs(FR[1]) < 1e-8 * W                           # roller carries no horizontal force
        # moment equilibrium about the origin (Y axis): reactions balance the line load
        # (resultant 2w at x = 3) and the self-weight (W at x = 3); the mesh is not
        # mirror-symmetric, so the bearings split the load only approximately evenly
        _, ML = total_reaction(model.results, brgL; about = zeros(3))
        _, MR = total_reaction(model.results, brgR; about = zeros(3))
        @test ML[2] + MR[2] ≈ -(3.0 * 2w + 3.0 * W) rtol = 1e-9
        @test FL[3] ≈ FR[3] rtol = 0.02
        # principal compression under the load (strut) and tension at mid-bottom (tie)
        g = ShellGeo(model)
        k_top = argmin(k -> norm(g.centers[:, k] - [3.0, 0, 2.7]), axes(g.centers, 2))
        k_bot = argmin(k -> norm(g.centers[:, k] - [3.0, 0, 0.2]), axes(g.centers, 2))
        @test g.N2[k_top] < 0
        @test g.N1[k_bot] > 0
    end
end
