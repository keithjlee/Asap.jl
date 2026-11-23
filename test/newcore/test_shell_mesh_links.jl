# Shell S5a: wall mesher, edge loads, LinkElement, and mixed wall–frame models.

using StaticArrays, LinearAlgebra
include("shell_test_utils.jl")

@testset "Wall mesher, edge loads, links (S5a)" begin
    conc = Material(3.0e7, 2.4, 0.2)
    wsec = ShellSection(conc, 0.2)

    @testset "mesher: topology, openings, placement" begin
        m = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 6.0, 3.0, wsec;
            openings = [(2.0, 0.0, 1.0, 2.1)], z_lines = [3.0], size = 0.5)
        @test all(op -> op isa ShellElement, m.elements)
        @test 2.0 in m.s && 3.0 in m.s && 2.1 in m.z
        # no element centre inside the door
        for el in m.elements
            c = sum(n -> n.position, el.nodes) / 4
            @test !(2.0 < c[1] < 3.0 && c[3] < 2.1)
        end
        @test sum(Asap.area, m.elements) ≈ 6.0 * 3.0 - 1.0 * 2.1
        # normal = along × up = −Y; e₁ along the wall
        R = local_frame(m.elements[1])
        @test R[3, :] ≈ [0, -1, 0] && R[1, :] ≈ [1, 0, 0]
        # bottom row supported except under the door
        @test all(n -> n.fixity == falses(6), filter(n -> n.position[3] == 0, m.nodes))
        @test_throws ArgumentError mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 6.0, 3.0, wsec;
            openings = [(5.5, 0.0, 1.0, 2.0)])
    end

    @testset "node reuse: L-shaped core shares its corner line" begin
        w1 = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 4.0, 3.0, wsec; size = 1.0)
        w2 = mesh_wall([4.0, 0, 0], [0.0, 1, 0], [0.0, 0, 1], 3.0, 3.0, wsec; size = 1.0, reuse = w1.nodes)
        shared = [n for n in w2.nodes if any(m -> m === n, w1.nodes)]
        @test length(shared) == length(w1.z)                          # the whole corner line
        nodes = unique(vcat(w1.nodes, w2.nodes))
        @test length(nodes) == length(w1.nodes) + length(w2.nodes) - length(shared)
        top = filter(n -> n.position[3] == 3.0, nodes)
        model = Model(nodes, AbstractElement{Float64}[w1.elements; w2.elements],
            AbstractLoad{Float64}[NodeForce(n, [10.0, 10.0, 0.0] ./ length(top)) for n in top])
        solve!(model)
        R = sum(n -> reaction(model.results, n)[SVector(1, 2, 3)], filter(n -> n.position[3] == 0, nodes))
        @test R ≈ [-10.0, -10.0, 0.0] atol = 1e-8
    end

    @testset "pier forces: storey shear splits between piers, exactly" begin
        m = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 6.0, 3.0, wsec;
            openings = [(2.0, 0.0, 1.0, 2.1)], z_lines = [1.0], size = 0.5)
        top = [m.grid[i, end] for i in eachindex(m.s)]
        model = Model(m.nodes, AbstractElement{Float64}[m.elements...],
            AbstractLoad{Float64}[NodeForce(n, [60.0 / length(top), 0, 0]) for n in top])
        solve!(model)
        eL, cL = wall_cut(m, 1.0; s_range = (0.0, 2.0))
        eR, cR = wall_cut(m, 1.0; s_range = (3.0, 6.0))
        FL, ML = section_cut(model, eL, cL)
        FR, MR = section_cut(model, eR, cR)
        @test FL[1] + FR[1] ≈ 60.0 rtol = 1e-10
        @test FR[1] > FL[1] > 0                      # the longer pier is stiffer
        eA, cA = wall_cut(m, 1.0)                    # whole-wall cut = sum of piers
        FA, _ = section_cut(model, eA, cA)
        @test FA ≈ FL + FR rtol = 1e-10
        @test_throws ArgumentError wall_cut(m, 1.234)
    end

    @testset "edge_loads: consistent line traction (uniform tension exact)" begin
        sec = ShellSection(Material(1000.0, 1.0, 0.3), 0.5)
        xy(i, j) = (i * 1.0 + (0 < i < 4 && j == 1 ? 0.2 : 0.0), 1.0j + (i == 2 && j == 1 ? 0.15 : 0.0))
        grid, els = membrane_mesh(4, 2, xy, sec; fix = (i, j) -> i == 0 ? (true, j == 0, true) : nothing)
        # elements in the last column, edge 2 (node 2 → 3) lies on x = 4
        right = [el for el in els if all(n -> n.position[1] >= 3.0 - 1e-12, (el.nodes[2], el.nodes[3])) &&
                                     el.nodes[2].position[1] == 4.0]
        loads = AbstractLoad{Float64}[l for el in right for l in edge_loads(el, 2, [2.0, 0.0, 0.0])]
        model = Model(vec(grid), els, loads)
        solve!(model)
        for el in els
            r = shell_resultants(model, el)
            Rp = r.frame[1:2, 1:2]
            @test Rp' * [r.N[1] r.N[3]; r.N[3] r.N[2]] * Rp ≈ [2.0 0; 0 0] atol = 1e-9
        end
    end

    @testset "LinkElement: stiff joint ≡ continuous wall; released joint" begin
        function split_wall(k)
            lower = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 2.0, 1.5, wsec; size = 0.5)
            upper = mesh_wall([0.0, 0, 1.5], [1.0, 0, 0], [0.0, 0, 1], 2.0, 1.5, wsec; size = 0.5,
                fixity_bottom = :free)                            # duplicate nodes at the joint
            joint = [LinkElement(lower.grid[i, end], upper.grid[i, 1], k) for i in eachindex(lower.s)]
            topn = [upper.grid[i, end] for i in eachindex(upper.s)]
            model = Model(vcat(lower.nodes, upper.nodes),
                AbstractElement{Float64}[lower.elements; upper.elements; joint],
                AbstractLoad{Float64}[NodeForce(n, [10.0 / length(topn), 0, 0]) for n in topn])
            solve!(model)
            return displacement(model.results, topn[end])[1]
        end
        cont = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 2.0, 3.0, wsec; size = 0.5)
        topc = [cont.grid[i, end] for i in eachindex(cont.s)]
        mc = Model(cont.nodes, AbstractElement{Float64}[cont.elements...],
            AbstractLoad{Float64}[NodeForce(n, [10.0 / length(topc), 0, 0]) for n in topc])
        solve!(mc)
        δc = displacement(mc.results, topc[end])[1]
        @test split_wall(fill(1e14, 6)) ≈ δc rtol = 1e-6
        # flexible horizontal joint (shear springs): more drift than the monolithic wall
        @test split_wall([1e6, 1e14, 1e14, 1e14, 1e14, 1e14]) > 1.01δc
        # zero-length enforced
        n1 = Node([0.0, 0, 0], :fixed); n2 = Node([0.0, 0, 1.0], :free)
        bad = LinkElement(n1, n2, fill(1.0, 6))
        @test_throws ArgumentError process!(Model([n1, n2], AbstractElement{Float64}[bad],
            AbstractLoad{Float64}[NodeForce(n2, [1.0, 0, 0])]))
        @test_throws ArgumentError LinkElement(n1, n2, zeros(6))
        @test !isempty(sprint(show, MIME"text/plain"(), bad))
        @test !isempty(sprint(show, MIME"text/plain"(), mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 1.0, 1.0, wsec)))
    end

    @testset "mixed model: in-plane beam framing into a wall transfers moment" begin
        # a cantilever beam framing into a thick (near-rigid) wall at a single
        # node, in the wall's plane: only the drilling DOF can hold the moment —
        # the tip deflection must approach the fixed-base cantilever PL³/3EI
        stiff = ShellSection(Material(3.0e10, 2.4, 0.2), 1.0)
        w = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 2.0, 4.0, stiff; z_lines = [2.0], size = 0.5)
        root = w.grid[end, findfirst(==(2.0), w.z)]                    # right edge, mid-height
        tip = Node([5.0, 0.0, 2.0], :free)
        bsec = Section(Material(2.0e8, 7.85, 0.3), 1e-2, 1e-4, 1e-4, 1e-5)
        beam = FrameElement(root, tip, bsec, :beam)
        P = 1.0
        model = Model(vcat(w.nodes, tip), AbstractElement{Float64}[w.elements..., beam],
            AbstractLoad{Float64}[NodeForce(tip, [0.0, 0.0, -P])])
        solve!(model)
        L = 3.0
        EI = 2.0e8 * 1e-4
        @test -displacement(model.results, tip)[3] ≈ P * L^3 / (3EI) rtol = 0.02
        # the beam's root moment equals the statics value P·L
        f = element_forces(model.results, beam)
        @test maximum(abs, f[SVector(4, 5, 6)]) ≈ P * L rtol = 1e-8
    end
end
