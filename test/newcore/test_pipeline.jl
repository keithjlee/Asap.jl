# Full-pipeline parity: rebuild every characterization model with the NEW
# core API, solve, and compare u / reactions / element forces / FEFs against
# the pinned legacy (v0.2.x) numerics. This is the Phase 1 exit criterion.
#
# DOF-space mapping: legacy truss models use 3 DOFs/node, the new core always
# 6 (rotations go inactive). Legacy element forces are 2-vectors for truss
# (axial at slots 1/2) vs the new uniform 12-vector (slots 1/7).

# new-core versions of the characterization model builders ------------------

function nc_truss_material(E)
    AsapNext.Material(E, 1.0, 1.0, 0.3)   # G/ρ/ν never consumed by truss kernels
end

function nc_truss_2d()
    N = AsapNext
    sec = N.Section(nc_truss_material(70.0), 4e3 / 1e6)
    # legacy TrussNode fixities, extended with free rotations (→ inactive)
    n1 = N.Node([0.0, 0.0, 0.0], [false, false, false, true, true, true])
    n2 = N.Node([10.0, 0.0, 0.0], [false, false, false, true, true, true])
    n3 = N.Node([0.0, 8.0, 0.0], [false, true, false, true, true, true])    # :yfree
    n4 = N.Node([6.0, 8.0, 0.0], [true, true, false, true, true, true])     # :free planarized later
    nodes = [n1, n2, n3, n4]
    elements = N.AbstractElement{Float64}[
        N.TrussElement(n1, n3, sec), N.TrussElement(n3, n4, sec),
        N.TrussElement(n1, n4, sec), N.TrussElement(n2, n3, sec),
        N.TrussElement(n2, n4, sec),
    ]
    loads = N.AbstractLoad{Float64}[
        N.NodeForce(n3, [0.0, -400.0, 0.0]),
        N.NodeForce(n4, [800.0, -400.0, 0.0]),
    ]
    return N.Model(nodes, elements, loads)
end

function nc_frame_2d()
    N = AsapNext
    mat = N.Material(29e3, 1.0, 1.0, 0.3)
    sec = N.Section(mat, 11.8, 310.0, 310.0, 1.0)
    n1 = N.Node([0.0, 0.0, 0.0], :fixed)
    n2 = N.Node([120.0, 240.0, 0.0], :free)
    n3 = N.Node([360.0, 240.0, 0.0], :fixed)
    e1 = N.FrameElement(n1, n2, sec; rollangle=0.0)
    e2 = N.FrameElement(n2, n3, sec; rollangle=0.0)
    model = N.Model([n1, n2, n3], N.AbstractElement{Float64}[e1, e2],
        N.AbstractLoad{Float64}[
            N.PointLoad(e1, 0.5, [0.0, -90.0, 0.0]),
            N.NodeMoment(n2, [0.0, 0.0, -125.0 * 12]),
            N.LineLoad(e2, [0.0, -1.5 / 12, 0.0]),
        ])
    N.planarize!(model)
    return model
end

function nc_truss_3d()
    N = AsapNext
    sec = N.Section(nc_truss_material(10e3), 8.4)
    freerot = [true, true, true]
    mk(pos, tfix) = N.Node(pos .* 12, vcat(tfix, freerot))
    n1 = mk([-6.0, 0.0, 8.0], [false, false, false])
    n2 = mk([12.0, 0.0, 8.0], [false, false, false])
    n3 = mk([6.0, 0.0, -8.0], [false, false, false])
    n4 = mk([-12.0, 0.0, -8.0], [false, false, false])
    n5 = mk([0.0, 24.0, 0.0], [true, true, true])
    els = N.AbstractElement{Float64}[N.TrussElement(a, n5, sec) for a in (n1, n2, n3, n4)]
    return N.Model([n1, n2, n3, n4, n5], els,
        N.AbstractLoad{Float64}[N.NodeForce(n5, [0.0, -100.0, -50.0])])
end

function nc_frame_3d()
    N = AsapNext
    mat = N.Material(29e3, 11.5e3, 1.0, 0.3)
    sec = N.Section(mat, 32.9, 716.0, 236.0, 15.1)
    n1 = N.Node([0.0, 0.0, 0.0], :free)
    n2 = N.Node([-240.0, 0.0, 0.0], :fixed)
    n3 = N.Node([0.0, -240.0, 0.0], :fixed)
    n4 = N.Node([0.0, 0.0, -240.0], :fixed)
    e1 = N.FrameElement(n2, n1, sec; rollangle=0.0)
    e2 = N.FrameElement(n3, n1, sec; rollangle=pi / 2)
    e3 = N.FrameElement(n4, n1, sec; rollangle=pi / 6)
    return N.Model([n1, n2, n3, n4], N.AbstractElement{Float64}[e1, e2, e3],
        N.AbstractLoad{Float64}[
            N.LineLoad(e1, [0.0, -3.0 / 12, 0.0]),
            N.NodeMoment(n1, [-150.0 * 12, 0.0, 150.0 * 12]),
        ])
end

function nc_cantilever(pointload::Bool)
    N = AsapNext
    mat = N.Material(30e6, 1.0, 1.0, 0.3)
    sec = N.Section(mat, 1.0, 57.1, 57.1, 1.0)
    n1 = N.Node([0.0, 0.0, 0.0], :fixed)
    n2 = N.Node([144.0, 0.0, 0.0], :free)
    e = N.FrameElement(n1, n2, sec)
    loads = pointload ?
            N.AbstractLoad{Float64}[N.PointLoad(e, 0.5, [0.0, -400.0, 0.0])] :
            N.AbstractLoad{Float64}[N.NodeForce(n2, [0.0, -400.0, 0.0])]
    model = N.Model([n1, n2], N.AbstractElement{Float64}[e], loads)
    N.planarize!(model)
    return model
end

function nc_portal_frame(release::Symbol)
    N = AsapNext
    H, L = 3000.0, 5000.0
    mat = N.Material(200.0, 77.0, 1.0, 0.3)
    sec = N.Section(mat, 1e4, 8e7, 3e7, 5e6)
    n1 = N.Node([0.0, 0.0, 0.0], :fixed)
    n2 = N.Node([0.0, 0.0, H], :free)
    n3 = N.Node([L, 0.0, H], :free)
    n4 = N.Node([L, 0.0, 0.0], :fixed)
    col1 = N.FrameElement(n1, n2, sec, :column)
    beam = N.FrameElement(n2, n3, sec, :beam; release=release)
    col2 = N.FrameElement(n4, n3, sec, :column)
    return N.Model([n1, n2, n3, n4], N.AbstractElement{Float64}[col1, beam, col2],
        N.AbstractLoad{Float64}[
            N.LineLoad(beam, [0.0, 0.0, -2.0]),
            N.PointLoad(beam, 0.4, [0.0, 3.0, -10.0]),
            N.NodeForce(n2, [5.0, 0.0, 0.0]),
        ])
end

# comparison helpers ---------------------------------------------------------

# expand a legacy per-node vector (3 or 6 DOFs/node) to the new 6-DOF space
function expand_legacy(v::Vector{Float64}, ndofs_legacy::Int, nnodes::Int)
    ndofs_legacy == 6 && return v
    out = zeros(6 * nnodes)
    for i in 1:nnodes
        out[6*(i-1).+(1:3)] = v[3*(i-1).+(1:3)]
    end
    return out
end

# legacy element force vector → new 12-vector layout
function expand_legacy_forces(f::Vector{Float64})
    length(f) == 12 && return f
    out = zeros(12)
    out[1] = f[1]
    out[7] = f[2]
    return out
end

function compare_solved(model, name::String, ndofs_legacy::Int)
    N = AsapNext
    N.solve!(model)
    res = model.results
    nn = length(model.nodes)

    @test res.u ≈ expand_legacy(FIXTURES["$name/u"], ndofs_legacy, nn) rtol = 1e-9 atol = 1e-11
    @test res.reactions ≈ expand_legacy(FIXTURES["$name/reactions"], ndofs_legacy, nn) rtol = 1e-9 atol = 1e-9
    @test res.compliance ≈ FIXTURES["$name/compliance"] rtol = 1e-9

    legacy_forces = FIXTURES["$name/element_forces"]
    for (i, el) in enumerate(model.elements)
        @test N.element_forces(res, el) ≈ expand_legacy_forces(legacy_forces[i]) rtol = 1e-8 atol = 1e-8
    end
end

# tests ----------------------------------------------------------------------

@testset "Pipeline parity vs pinned legacy numerics" begin
    N = AsapNext

    @testset "FEF engine ≡ pinned q_local/q: $r" for r in CHAR_RELEASES
        sec = fef_section()
        ends = N.EndConditions(r)
        x1, x2 = FEF_X1, FEF_X2
        L = N.element_length(x1, x2)
        rollangle = pi / 2

        dummy = N.Node(x1, :fixed)
        dummy2 = N.Node(x2, :fixed)
        el = N.FrameElement(dummy, dummy2, sec; release=r)

        line = N.LineLoad(el, [0.4, -1.2, 0.9])
        point = N.PointLoad(el, 0.3, [110.0, -70.0, 25.0])

        for (load, kind) in ((line, "lineload"), (point, "pointload"))
            q = N.fixed_end_forces(load, sec, ends, x1, x2, rollangle)
            q̃ = N.condense_fef(q, sec, L, ends)

            # DELIBERATE DEVIATION (documented in docs/MODERNIZATION.md):
            # legacy q_local(::PointLoad) splits the axial component 50/50
            # regardless of load position; the consistent (and statically
            # exact) distribution is by lever rule: −P·(1−ξ) / −P·ξ. All
            # non-axial slots must still match the pins exactly.
            nonaxial = setdiff(1:12, (1, 7))
            for (qq, key) in ((q, "q_local"), (q̃, "q"))
                @test Vector(qq)[nonaxial] ≈ FIXTURES["fef_$r/$key/$kind"][nonaxial] rtol = 1e-10 atol = 1e-10
            end
            if kind == "pointload"
                Λ = N.local_frame(x1, x2, rollangle)
                Pax = (Λ*load.value)[1]
                ξ = load.position
                @test q[1] ≈ -Pax * (1 - ξ) rtol = 1e-12
                @test q[7] ≈ -Pax * ξ rtol = 1e-12
                # totals always agree with legacy (statics preserved)
                @test q[1] + q[7] ≈ FIXTURES["fef_$r/q_local/$kind"][1] +
                                    FIXTURES["fef_$r/q_local/$kind"][7] rtol = 1e-10
            else
                @test Vector(q)[[1, 7]] ≈ FIXTURES["fef_$r/q_local/$kind"][[1, 7]] rtol = 1e-10 atol = 1e-10
            end
        end
    end

    @testset "solve parity: truss_2d" begin
        compare_solved(nc_truss_2d(), "truss_2d", 3)
    end
    @testset "solve parity: truss_3d" begin
        compare_solved(nc_truss_3d(), "truss_3d", 3)
    end
    @testset "solve parity: frame_2d" begin
        compare_solved(nc_frame_2d(), "frame_2d", 6)
    end
    @testset "solve parity: frame_3d" begin
        compare_solved(nc_frame_3d(), "frame_3d", 6)
    end
    @testset "solve parity: cantilevers" begin
        compare_solved(nc_cantilever(false), "cantilever_tip", 6)
        compare_solved(nc_cantilever(true), "cantilever_midpoint", 6)
    end
    @testset "solve parity: portal_$r" for r in CHAR_RELEASES
        compare_solved(nc_portal_frame(r), "portal_$r", 6)
    end

    @testset "truss-only model excludes rotations structurally" begin
        model = nc_truss_3d()
        N.solve!(model)
        part = model.cache.partition
        # every rotational DOF is inactive — the solve was translations-only
        @test length(part.inactive) == 3 * length(model.nodes)
        @test all(g -> mod1(g, 6) > 3, part.inactive)
    end

    @testset "mixed truss + frame in ONE model" begin
        # a frame cantilever with a truss tie-rod to ground — impossible in
        # the legacy library, the raison d'être of per-element DOF activity
        mat = N.Material(200.0, 77.0, 1.0, 0.3)
        fsec = N.Section(mat, 1e4, 8e7, 3e7, 5e6)
        tsec = N.Section(mat, 1e3)

        n1 = N.Node([0.0, 0.0, 0.0], :fixed)
        n2 = N.Node([3000.0, 0.0, 0.0], :free)
        n3 = N.Node([3000.0, 0.0, 3000.0], :fixed)   # tie anchor

        beam = N.FrameElement(n1, n2, fsec, :beam)
        tie = N.TrussElement(n2, n3, tsec, :tie)
        load = N.NodeForce(n2, [0.0, 0.0, -100.0])

        model = N.Model([n1, n2, n3], N.AbstractElement{Float64}[beam, tie],
            N.AbstractLoad{Float64}[load])
        N.solve!(model)
        res = model.results

        # the tie carries tension; the beam tip deflects far less than alone
        @test N.axial_force(res, tie) > 0
        u_with = N.displacement(res, n2)[3]

        model2 = N.Model([n1, n2, n3], N.AbstractElement{Float64}[
                N.FrameElement(n1, n2, fsec, :beam)],
            N.AbstractLoad{Float64}[N.NodeForce(n2, [0.0, 0.0, -100.0])])
        N.solve!(model2)
        u_alone = N.displacement(model2.results, n2)[3]
        @test abs(u_with) < abs(u_alone) / 2

        # equilibrium: vertical reactions balance the load
        Rz = sum(N.reaction(res, n)[3] for n in model.nodes)
        @test Rz ≈ 100.0 rtol = 1e-9
    end

    @testset "spring support" begin
        # cantilever with a vertical spring at the tip: u = P/(k_beam + k_spring)
        mat = N.Material(30e6, 1.0, 1.0, 0.3)
        sec = N.Section(mat, 1.0, 57.1, 57.1, 1.0)
        n1 = N.Node([0.0, 0.0, 0.0], :fixed)
        n2 = N.Node([144.0, 0.0, 0.0], :free)
        e = N.FrameElement(n1, n2, sec)
        P = 400.0
        kbeam = 3 * N.EIx(sec) / 144.0^3          # tip stiffness of a cantilever
        kspring = 2 * kbeam

        model = N.Model([n1, n2], N.AbstractElement{Float64}[e],
            N.AbstractLoad{Float64}[N.NodeForce(n2, [0.0, -P, 0.0])];
            springs=[N.NodalSpring(n2, [0.0, kspring, 0.0, 0.0, 0.0, 0.0])])
        N.planarize!(model)
        N.solve!(model)

        u = N.displacement(model.results, n2)[2]
        @test u ≈ -P / (kbeam + kspring) rtol = 1e-6
    end

    @testset "inactive-DOF load errors are diagnostic" begin
        model = nc_truss_2d()
        push!(model.loads, N.NodeMoment(model.nodes[4], [0.0, 0.0, 100.0]))
        err = try
            N.solve!(model; reprocess=true)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("inactive", err.msg) && occursin("truss", err.msg)
    end
end

@testset "zero-rigidity guard (process!-time validation)" begin
    N = AsapNext
    mat = N.Material(200e6, 77e6, 8.0, 0.3)
    axial_only = N.Section(mat, 1e-3)                 # Ix = Iy = J = 0
    full = N.Section(mat, 1e-3, 1e-5, 1e-5, 1e-6)
    n1 = N.Node([0.0, 0.0, 0.0], :fixed)
    n2 = N.Node([2.0, 0.0, 0.0], :fixed)

    # axial-only section on a rigid-ended frame element: caught with a clear error
    bad = N.FrameElement(n1, n2, axial_only)
    m = N.Model([n1, n2], N.AbstractElement{Float64}[bad],
        N.AbstractLoad{Float64}[N.NodeForce(n2, [0.0, 0.0, -1.0])])
    @test_throws ArgumentError N.process!(m)

    # ...but fully released (:freefree brace pattern) is legitimate
    a1 = N.Node([0.0, 0.0, 0.0], :fixed); a2 = N.Node([2.0, 0.0, 0.0], :fixed)
    brace = N.FrameElement(a1, a2, axial_only; release = :freefree)
    m2 = N.Model([a1, a2], N.AbstractElement{Float64}[brace],
        N.AbstractLoad{Float64}[N.NodeForce(a2, [0.0, 0.0, -1.0])])
    @test N.process!(m2) isa N.Model

    # :joist keeps torsion engaged, so J = 0 is caught
    b1 = N.Node([0.0, 0.0, 0.0], :fixed); b2 = N.Node([2.0, 0.0, 0.0], :fixed)
    joist_noJ = N.FrameElement(b1, b2, N.Section(mat, 1e-3, 1e-5, 1e-5, 0.0);
        release = :joist)
    m3 = N.Model([b1, b2], N.AbstractElement{Float64}[joist_noJ],
        N.AbstractLoad{Float64}[N.NodeForce(b2, [0.0, 0.0, -1.0])])
    @test_throws ArgumentError N.process!(m3)

    # truss with EA = 0 is caught; with EA > 0 passes
    c1 = N.Node([0.0, 0.0, 0.0], :pinned); c2 = N.Node([2.0, 0.0, 0.0], :pinned)
    dead_bar = N.TrussElement(c1, c2, N.Section(mat, 0.0))
    m4 = N.Model([c1, c2], N.AbstractElement{Float64}[dead_bar],
        N.AbstractLoad{Float64}[N.NodeForce(c2, [1.0, 0.0, 0.0])])
    @test_throws ArgumentError N.process!(m4)

    # VariableElement: every segment needs full rigidities
    d1 = N.Node([0.0, 0.0, 0.0], :fixed); d2 = N.Node([4.0, 0.0, 0.0], :fixed)
    varbad = N.VariableElement(d1, d2,
        N.AbstractSection{Float64}[full, axial_only], [0.5])
    m5 = N.Model([d1, d2], N.AbstractElement{Float64}[varbad],
        N.AbstractLoad{Float64}[N.NodeForce(d2, [0.0, 0.0, -1.0])])
    @test_throws ArgumentError N.process!(m5)
end

@testset "assembly is allocation-free (truss AND frame groups)" begin
    N = AsapNext
    sec = N.Section(N.Material(200e6, 77e6, 8.0, 0.3), 1e-3, 1e-5, 1e-5, 1e-6)

    #truss group: regression for the hvcat heap fallback in truss_stiffness
    #(this once allocated MBs per assembly and was 8× slower on Julia 1.12)
    #tetrahedron: three pinned base nodes + loaded apex — stable in 3D
    tn = [N.Node([0.0, 0.0, 0.0], :pinned), N.Node([1.0, 0.0, 0.0], :pinned),
        N.Node([0.5, 1.0, 0.0], :pinned), N.Node([0.5, 0.4, 1.0], :free)]
    tels = N.AbstractElement{Float64}[N.TrussElement(tn[i], tn[4], sec) for i in 1:3]
    tm = N.Model(tn, tels, N.AbstractLoad{Float64}[N.NodeForce(tn[4], [0.0, -1.0, -1.0])])
    N.solve!(tm)
    N.assemble_K!(tm.cache, tm)                       # warm up
    @test @allocated(N.assemble_K!(tm.cache, tm)) <= 256

    #frame group
    fn = [N.Node([Float64(i), 0.0, 0.0], i == 1 ? :fixed : :free) for i in 1:3]
    fels = N.AbstractElement{Float64}[N.FrameElement(fn[i], fn[i+1], sec) for i in 1:2]
    fm = N.Model(fn, fels, N.AbstractLoad{Float64}[N.NodeForce(fn[3], [0.0, 0.0, -1.0])])
    N.solve!(fm)
    N.assemble_K!(fm.cache, fm)
    @test @allocated(N.assemble_K!(fm.cache, fm)) <= 256
end

@testset "mutate-and-resolve freshness + factorization reuse" begin
    N = AsapNext
    mat = N.Material(200e6, 77e6, 8.0, 0.3)
    mk(k, A) = begin
        n1 = N.Node([0.0, 0.0, 0.0], :fixed)
        n2 = N.Node([0.0, 0.0, 3.0], :free)
        el = N.FrameElement(n1, n2, N.Section(mat, A, 1e-4, 5e-5, 1e-6))
        sp = N.NodalSpring(n2, [k, 0.0, 0.0, 0.0, 0.0, 0.0])
        N.Model([n1, n2], N.AbstractElement{Float64}[el],
            N.AbstractLoad{Float64}[N.NodeForce(n2, [10.0, 0.0, 0.0])]; springs = [sp])
    end

    m = mk(1e4, 1e-2)
    N.solve!(m)
    fc = m.cache.factorization
    @test fc isa N.FactorizationCache

    # replace spring (same pattern, new VALUE) + section; resolve WITHOUT reprocess
    m.springs[1] = N.NodalSpring(m.nodes[2], [5e4, 0.0, 0.0, 0.0, 0.0, 0.0])
    m.elements[1].section = N.Section(mat, 2e-2, 1e-4, 5e-5, 1e-6)
    N.solve!(m)
    fresh = mk(5e4, 2e-2)
    N.solve!(fresh)
    @test N.displacement(m.results, m.nodes[2]) ≈ N.displacement(fresh.results, fresh.nodes[2]) rtol = 1e-12

    # the factorization object was REUSED (numeric-only refactorization)
    @test m.cache.factorization === fc

    # pure path reads replaced spring values too
    st = N.extract_state(m)
    K = N.assemble_K(m.cache, st)
    @test K ≈ m.cache.K rtol = 1e-12

    # sticky solver: default resolves keep the same FactorizationCache
    N.solve!(m)
    @test m.cache.factorization === fc
end
