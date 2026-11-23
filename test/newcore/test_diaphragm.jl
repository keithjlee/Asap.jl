# Shell S5b: rigid diaphragms (exact transformation u = C·û).
# Oracle: closed-form storey stiffness of a rigid floor on cantilever columns
# (lateral 3EI/H³ each way, torsion GJ/H) — exact for Euler–Bernoulli frames.

using StaticArrays, LinearAlgebra, Zygote, FiniteDifferences

function _storey(; H = 3.0, Is = (1.0e-4, 2.0e-4, 1.5e-4, 3.0e-4), load = (10.0, -4.0, 25.0))
    mat = Material(2.0e8, 7.85, 0.3)
    xy = ((0.0, 0.0), (6.0, 0.0), (0.0, 4.0), (6.0, 4.0))
    base = [Node([x, y, 0.0], :fixed) for (x, y) in xy]
    top = [Node([x, y, H], :free) for (x, y) in xy]
    secs = [Section(mat, 1e-2, I, I, 0.5I) for I in Is]
    cols = AbstractElement{Float64}[FrameElement(base[k], top[k], secs[k], :col) for k in 1:4]
    d = RigidDiaphragm(top; id = :roof)
    Fx, Fy, Mz = load
    m = Model(vcat(base, top), cols, AbstractLoad{Float64}[NodeForce(d.master, [Fx, Fy, 0.0]),
        NodeMoment(d.master, [0.0, 0.0, Mz])]; diaphragms = [d])
    return m, d, top, secs, mat, H
end

@testset "Rigid diaphragms (S5b)" begin

    @testset "storey torsion: exact against the 3×3 rigid-floor stiffness" begin
        m, d, top, secs, mat, H = _storey()
        solve!(m)
        xm = d.master.position
        K = zeros(3, 3)
        for (k, n) in enumerate(top)
            ki = 3 * EIx(secs[k]) / H^3
            dx, dy = n.position[1] - xm[1], n.position[2] - xm[2]
            K .+= ki .* [1 0 -dy; 0 1 dx; -dy dx dx^2+dy^2]
            K[3, 3] += GJ(secs[k]) / H
        end
        uref = K \ [10.0, -4.0, 25.0]
        um = displacement(m.results, d.master)
        @test [um[1], um[2], um[6]] ≈ uref rtol = 1e-9
        # slaves follow the rigid-body field exactly
        for n in top
            us = displacement(m.results, n)
            r = n.position - xm
            @test us[1] ≈ um[1] - um[6] * r[2] rtol = 1e-9
            @test us[2] ≈ um[2] + um[6] * r[1] rtol = 1e-9
            @test us[6] ≈ um[6] rtol = 1e-9
        end
        # equilibrium: base reactions balance the applied storey forces and torque
        R = sum(n -> reaction(m.results, n), filter(n -> n.position[3] == 0, m.nodes))
        @test R[SVector(1, 2)] ≈ [-10.0, 4.0] rtol = 1e-9
        @test m.results.compliance ≈ dot([10.0, -4.0, 25.0], uref) rtol = 1e-9
        # the auto master: in-plane free, out-of-plane fixed
        @test d.master.fixity == Bool[1, 1, 0, 0, 0, 1]
        @test any(n -> n === d.master, m.nodes)
    end

    @testset "symmetric shear walls share the storey shear equally" begin
        sec = ShellSection(Material(3.0e7, 2.4, 0.2), 0.2)
        w1 = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 4.0, 3.0, sec; size = 1.0)
        w2 = mesh_wall([0.0, 6.0, 0], [1.0, 0, 0], [0.0, 0, 1], 4.0, 3.0, sec; size = 1.0)
        topn = vcat([w1.grid[i, end] for i in eachindex(w1.s)], [w2.grid[i, end] for i in eachindex(w2.s)])
        d = RigidDiaphragm(topn)
        for (Fx, Mz) in ((100.0, 0.0), (100.0, 60.0))
            m = Model(vcat(w1.nodes, w2.nodes), AbstractElement{Float64}[w1.elements; w2.elements],
                AbstractLoad{Float64}[NodeForce(d.master, [Fx, 0, 0]), NodeMoment(d.master, [0, 0, Mz])];
                diaphragms = [d])
            solve!(m)
            V1 = -sum(reaction(m.results, n)[1] for n in w1.nodes if n.position[3] == 0)
            V2 = -sum(reaction(m.results, n)[1] for n in w2.nodes if n.position[3] == 0)
            @test V1 + V2 ≈ Fx rtol = 1e-9
            # torque: mostly the in-plane wall-shear couple (Mz ≈ 3(V₁ − V₂)); the
            # walls' out-of-plane bending/torsion carries the remaining few tenths of
            # a percent, so check the EXACT global balance of all base reactions
            @test (V1 - V2) * 3.0 ≈ Mz atol = 0.01 * abs(Mz) + 1e-9 * Fx
            xm = d.master.position
            Mr = sum(n -> cross(n.position - xm, reaction(m.results, n)[SVector(1, 2, 3)])[3] +
                          reaction(m.results, n)[6], filter(n -> n.position[3] == 0, m.nodes))
            @test Mr + Mz ≈ 0 atol = 1e-8 * Fx
        end
    end

    @testset "cases, pure-path parity and gradient (offsets frozen)" begin
        m, d, top, secs, mat, H = _storey()
        push!(m.loads, NodeForce(d.master, [0.0, 7.0, 0.0]; case = :W))
        cr = solve_cases!(m; cases = [:LC1, :W])
        solve!(m)
        @test cr.results[1].u + cr.results[2].u ≈ m.results.u rtol = 1e-10
        st = extract_state(m)
        @test compliance(m, st) ≈ m.results.compliance rtol = 1e-10
        @test solve(m, st) ≈ m.results.u rtol = 1e-10
        # gradient w.r.t. column-top heights (the master's offsets stay frozen)
        f(X) = compliance(m, ModelState{Float64}(X, st.sections))
        gz = Zygote.gradient(f, st.X)[1]
        gfd = FiniteDifferences.grad(central_fdm(5, 1), f, st.X)[1]
        @test gz ≈ gfd rtol = 1e-5
    end

    @testset "validation" begin
        n1 = Node([0.0, 0, 3], :pinned)                  # in-plane supported slave
        n2 = Node([2.0, 0, 3], :free)
        b1 = Node([0.0, 0, 0], :fixed); b2 = Node([2.0, 0, 0], :fixed)
        sec = Section(Material(2.0e8, 7.85, 0.3), 1e-2, 1e-4, 1e-4, 1e-5)
        els = AbstractElement{Float64}[FrameElement(b1, n1, sec), FrameElement(b2, n2, sec)]
        d = RigidDiaphragm([n1, n2])
        @test_throws ArgumentError process!(Model([b1, b2, n1, n2], els,
            AbstractLoad{Float64}[NodeForce(d.master, [1.0, 0, 0])]; diaphragms = [d]))
        @test_throws ArgumentError RigidDiaphragm([n1])
        @test_throws ArgumentError RigidDiaphragm([n1, n2]; axis = :W)
        @test_throws ArgumentError RigidDiaphragm([n1, n2]; master = n1)
        @test !isempty(sprint(show, MIME"text/plain"(), d)) && !isempty(sprint(show, d))
    end
end
