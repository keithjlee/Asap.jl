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
        # holes may touch the outline or each other (vertex, vertex on an edge, shared
        # edge: notches drawn as holes, adjacent openings as two); the mesh covers
        # exactly the material, with the same quality
        sq = [(0.0, 0.0), (6.0, 0.0), (6.0, 3.0), (0.0, 3.0)]
        box(x0, y0, x1, y1) = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        for (holes_, A) in (
                ([box(2.0, 0.0, 3.0, 1.0)], 17.0),                                     # notch
                ([[(0.0, 0.0), (1.0, 0.5), (0.5, 1.0)]], 17.625),                      # at the outline's corner
                ([[(3.0, 0.0), (3.5, 0.8), (2.5, 0.8)]], 17.6),                        # vertex on the outline's edge
                ([[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0)], [(2.0, 2.0), (3.0, 2.0), (3.0, 2.5)]], 17.25),   # holes share a vertex
                ([box(1.0, 1.0, 2.0, 2.0), box(2.0, 1.0, 3.0, 2.0)], 16.0),            # holes share an edge
                ([box(1.0, 1.0, 2.0, 2.0), box(2.0, 1.5, 3.0, 2.5)], 16.0))            # part of an edge
            mt = mesh_polygon(sq, sec; holes = holes_, h = 0.25)
            @test sum(Asap.area, mt.elements) ≈ A rtol = 1e-10
            @test minimum(mt.elements) do el
                x = map(n -> n.position, el.nodes)
                minimum(k -> acosd(clamp(dot(normalize(x[mod1(k + 1, 3)] - x[k]), normalize(x[mod1(k + 2, 3)] - x[k])), -1, 1)), 1:3)
            end >= 25 - 1e-6
        end
        # a non-convex (L-shaped) outline with a hole touching its re-entrant corner
        Lsh = [(0.0, 0.0), (4.0, 0.0), (4.0, 1.0), (1.0, 1.0), (1.0, 3.0), (0.0, 3.0)]
        mL = mesh_polygon(Lsh, sec; holes = [box(1.0, 0.6, 1.6, 1.0)], h = 0.2)
        @test sum(Asap.area, mL.elements) ≈ 6.0 - 0.24 rtol = 1e-10
        # a notch's mouth lies on the outline: no material on either side of it, so only
        # its two end nodes are mesh nodes (holes are stored clockwise: edge order reversed)
        mt = mesh_polygon(sq, sec; holes = [box(2.0, 0.0, 3.0, 1.0)], h = 0.25)
        mouth = only(filter(k -> all(n -> n.position[2] ≈ 0, boundary_nodes(mt, 2, k)), 1:4))
        @test [n.position[1] for n in boundary_nodes(mt, 2, mouth)] ≈ [3.0, 2.0]
        # refused: a hole crossing the outline, overlapping holes, a hole outside
        @test_throws ArgumentError mesh_polygon(sq, sec; holes = [box(5.0, 1.0, 7.0, 2.0)], h = 0.5)
        @test_throws ArgumentError mesh_polygon(sq, sec; holes = [box(1.0, 1.0, 2.0, 2.0), box(1.5, 1.5, 2.5, 2.5)], h = 0.5)
        @test_throws ArgumentError mesh_polygon(sq, sec; holes = [box(7.0, 1.0, 8.0, 2.0)], h = 0.5)
        # reproducible: same input and seed → identical mesh
        m2 = mesh_polygon(outline, sec; holes = [circle, L], h = 0.25, min_angle = 25)
        @test m2.points == m.points
        # an extra point becomes a node
        mp = mesh_polygon(outline, sec; h = 0.5, points = [(3.0, 1.5)])
        @test any(n -> n.position ≈ [3.0, 1.5, 0.0], mp.nodes)
    end

    @testset "read_gmsh: MSH 4.1 ASCII, physical groups, planar and curved" begin
        # 2 × 1 panel: one quad (1,2,5,6) and two triangles; physical curve "bottom"
        # (y = 0) and physical surface "panel"; nodes shifted by `shift(x, y)`
        function msh(z)
            P = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (0.0, 1.0)]
            nodes = join(["$(x) $(y) $(z(x, y))" for (x, y) in P], "\n")
            """
            \$MeshFormat
            4.1 0 8
            \$EndMeshFormat
            \$PhysicalNames
            2
            1 1 "bottom"
            2 2 "panel"
            \$EndPhysicalNames
            \$Entities
            0 1 1 0
            1 0 0 0 2 0 0 1 1 0
            1 0 0 0 2 1 1 1 2 0
            \$EndEntities
            \$Nodes
            1 6 1 6
            2 1 0 6
            1
            2
            3
            4
            5
            6
            $nodes
            \$EndNodes
            \$Elements
            3 6 1 6
            1 1 1 2
            1 1 2
            2 2 3
            2 1 3 1
            3 1 2 5 6
            2 1 2 2
            4 2 3 4
            5 2 4 5
            \$EndElements
            """
        end
        f = joinpath(mktempdir(), "panel.msh")
        write(f, msh((x, y) -> 0.0))
        m = read_gmsh(f, sec)
        @test size(m.points) == (2, 6)
        @test count(el -> el isa ShellElement{Float64,4}, m.elements) == 1
        @test count(el -> el isa ShellElement{Float64,3}, m.elements) == 2
        @test sort([n.position[1] for n in group_nodes(m, "bottom")]) ≈ [0.0, 1.0, 2.0]
        @test length(group_elements(m, "panel")) == 3
        @test sum(Asap.area, m.elements) ≈ 2.0
        @test all(el -> local_frame(el)[3, :] ≈ local_frame(m.elements[1])[3, :], m.elements)
        @test_throws ArgumentError group_nodes(m, "top")
        # the same panel tilted out of the XY plane: still planar, lengths preserved
        write(f, msh((x, y) -> 0.5y))
        mt = read_gmsh(f, sec)
        @test size(mt.points) == (2, 6)
        @test sum(Asap.area, mt.elements) ≈ 2.0 * sqrt(1.25)
        @test length(nodes_on_segment(mt, mt.points[:, 1], mt.points[:, 3])) == 3
        # curved (z = x²/4): a spatial mesh with 3D coordinates; segment queries in 3D
        write(f, msh((x, y) -> x^2 / 4))
        mc = read_gmsh(f, sec)
        @test size(mc.points) == (3, 6)
        @test length(nodes_on_segment(mc, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0])) == 2
        @test_throws ArgumentError read_gmsh(f, sec; planar = true)
        @test_throws ArgumentError elements_in(mc, [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])
        # solve on the curved mesh: bottom edge clamped, a load at the far corner
        foreach(n -> fixnode!(n, :fixed), group_nodes(mc, "bottom"))
        corner = only(filter(n -> n.position[1] ≈ 2 && n.position[2] ≈ 1, mc.nodes))
        model = Model(mc.nodes, AbstractElement{Float64}[mc.elements...],
            AbstractLoad{Float64}[NodeForce(corner, [0.0, 0.0, -1.0])])
        solve!(model)
        F, _ = total_reaction(model.results, group_nodes(mc, "bottom"))
        @test F ≈ [0.0, 0.0, 1.0] atol = 1e-9
        # unsupported formats are refused clearly
        write(f, replace(msh((x, y) -> 0.0), "4.1 0 8" => "2.2 0 8"))
        @test_throws ArgumentError read_gmsh(f, sec)
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

    @testset "node reuse: a triangulated wall joined to a frame" begin
        # a 6 × 3 wall in the XZ plane; a beam framing into its top edge at x = 2 and
        # x = 4.5, and a column joint inside the panel at (3, 1.5): all become mesh
        # vertices and the mesh uses those very Node objects
        steel = Section(Material(2.0e8, 7.7e7, 7.85, 0.3), 1e-2, 1e-4, 1e-4, 1e-5)
        a = Node([2.0, 0.0, 3.0], :free)
        b = Node([4.5, 0.0, 3.0], :free)
        c = Node([3.0, 0.0, 1.5], :free)
        off = Node([3.0, 2.0, 1.5], :free)                 # off the panel's plane: ignored
        tip = Node([4.5, 3.0, 3.0], :free)
        outline = [(0.0, 0.0), (6.0, 0.0), (6.0, 3.0), (0.0, 3.0)]
        m = mesh_polygon(outline, sec; h = 0.4, e1 = [1.0, 0, 0], e2 = [0.0, 0, 1], reuse = [a, b, c, off])
        @test any(n -> n === a, m.nodes) && any(n -> n === b, m.nodes) && any(n -> n === c, m.nodes)
        @test !any(n -> n === off, m.nodes)
        foreach(n -> fixnode!(n, :fixed), nodes_on_segment(m, (0.0, 0.0), (6.0, 0.0)))
        frame = AbstractElement{Float64}[FrameElement(b, tip, steel), FrameElement(c, off, steel)]
        nodes = unique([a, b, c, off, tip, m.nodes...])
        loads = AbstractLoad{Float64}[NodeForce(tip, [0.0, 0.0, -10.0]), NodeForce(off, [5.0, 0.0, 0.0])]
        model = Model(nodes, AbstractElement{Float64}[m.elements..., frame...], loads)
        solve!(model)
        F, _ = total_reaction(model.results, nodes_on_segment(m, (0.0, 0.0), (6.0, 0.0)))
        @test F ≈ [-5.0, 0.0, 10.0] rtol = 1e-9                # the frame loads reach the wall's base
        ub = displacement(model.results, b)
        # listing a reused node twice is harmless (as with mesh_wall's reuse): its first
        # slot connects to nothing and is retired as inactive; the results are the same
        dup = Model([a, b, c, off, tip, m.nodes...], AbstractElement{Float64}[m.elements..., frame...], loads)
        solve!(dup)
        @test displacement(dup.results, b) ≈ ub rtol = 1e-10
        # array constructor: reuse by position
        P = [0.0 1.0 1.0 0.0; 0.0 0.0 1.0 1.0]
        corner = Node([1.0, 1.0, 0.0], :fixed)
        ma = shell_mesh(P, [(1, 2, 3), (1, 3, 4)], sec; reuse = [corner])
        @test any(n -> n === corner, ma.nodes)
        @test count(n -> n.position ≈ [1.0, 1.0, 0.0], ma.nodes) == 1
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
