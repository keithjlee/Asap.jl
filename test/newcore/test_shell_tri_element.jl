# S7d: the triangular ShellElement in assembled models — convergence on the standard
# membrane, plate and shell benchmarks (meshes from the quad utilities, each quad split
# into two triangles), patch tests on irregular meshes, mixed quad/triangle meshes,
# loads, resultants, equilibrium, plotting arrays and pure-path parity.

using LinearAlgebra, StaticArrays
include(joinpath(@__DIR__, "shell_test_utils.jl"))

# split each quad (1,2,3,4) into (1,2,3) + (1,3,4), or along 2–4 when `alt`
_tri_split(els; alt = false) = AbstractElement{Float64}[
    t for el in els for t in (alt ?
        (ShellElement((el.nodes[1], el.nodes[2], el.nodes[4]), el.section),
         ShellElement((el.nodes[2], el.nodes[3], el.nodes[4]), el.section)) :
        (ShellElement((el.nodes[1], el.nodes[2], el.nodes[3]), el.section),
         ShellElement((el.nodes[1], el.nodes[3], el.nodes[4]), el.section)))]

function _tri_lumped(els, q)
    acc = Dict{Node{Float64},Vector{Float64}}()
    for el in els, n in el.nodes
        acc[n] = get(acc, n, zeros(3)) .+ q .* (Asap.area(el) / length(el.nodes))
    end
    return AbstractLoad{Float64}[NodeForce(n, f) for (n, f) in acc]
end

function _cook_tri(n)
    sec = ShellSection(Material(1.0, 1.0, 1 / 3), 1.0)
    grid, els = membrane_mesh(n, n, (i, j) -> (48i / n, 44i / n + (44 - 28i / n) * j / n), sec;
        fix = (i, j) -> i == 0 ? (true, true, true) : nothing)
    loads = [NodeForce(grid[n+1, j+1], [0.0, (j == 0 || j == n ? 0.5 : 1.0) / n, 0.0]) for j in 0:n]
    m = Model(vec(grid), _tri_split(els), AbstractLoad{Float64}[loads...])
    solve!(m)
    return displacement(m.results, grid[n+1, n÷2+1])[2]
end

function _plate_tri(n; edge = :ss, a = 1.0, t = 0.01, E = 1.0e4, ν = 0.3)
    sec = ShellSection(Material(E, 1.0, ν), t)
    onedge(i, j) = i == 0 || j == 0 || i == n || j == n
    fix(i, j) = onedge(i, j) ? (edge == :ss ? Bool[0, 0, 0, 1, 1, 0] : Bool[0, 0, 0, 0, 0, 0]) :
                Bool[0, 0, 1, 1, 1, 0]
    grid, els = grid_shell_mesh(n, n, (i, j) -> (a * i / n, a * j / n, 0.0), sec; fix = fix)
    tris = _tri_split(els)
    m = Model(vec(grid), tris, _tri_lumped(tris, [0.0, 0.0, -1.0]))
    solve!(m)
    D = E * t^3 / (12 * (1 - ν^2))
    return -displacement(m.results, grid[n÷2+1, n÷2+1])[3] * D / a^4
end

function _scordelis_tri(n)
    R, L, Φ = 25.0, 50.0, deg2rad(40.0)
    sec = ShellSection(Material(4.32e8, 1.0, 0.0), 0.25)
    X(i, j) = (L / 2 * i / n, R * sin(Φ * j / n), R * cos(Φ * j / n))
    function fix(i, j)
        f = trues(6)
        i == 0 && (f[2] = f[3] = false)
        i == n && (f[1] = f[5] = f[6] = false)
        j == 0 && (f[2] = f[4] = f[6] = false)
        return collect(f)
    end
    grid, els = grid_shell_mesh(n, n, X, sec; fix = fix)
    tris = _tri_split(els)
    m = Model(vec(grid), tris, _tri_lumped(tris, [0.0, 0.0, -90.0]))
    solve!(m)
    return -displacement(m.results, grid[n+1, n+1])[3]
end

@testset "Triangular ShellElement (S7d)" begin

    @testset "construction, frame, validation" begin
        sec = ShellSection(Material(3.0e4, 2.4e-3, 0.2), 0.2)
        n1, n2, n3 = Node([0.0, 0, 0], :free), Node([2.0, 0, 0.5], :free), Node([0.3, 1.5, 0.2], :free)
        el = ShellElement((n1, n2, n3), sec)
        @test el isa ShellElement{Float64,3}
        @test Asap.ndofs(el) == 18
        R = local_frame(el)
        @test R * R' ≈ I(3) atol = 1e-14
        @test R[1, :] ≈ normalize([2.0, 0, 0.5])
        @test Asap.area(el) ≈ norm(cross([2.0, 0, 0.5], [0.3, 1.5, 0.2])) / 2
        bad = ShellElement((n1, n2, Node([4.0, 0, 1.0], :free)), sec)        # collinear
        @test_throws ArgumentError Asap.validate_rigidities(bad)
        @test_throws ArgumentError ShellElement((n1, n2), sec)
        # global kernel: rigid-body null space of rank 6 in any orientation
        K = shell_stiffness(sec, n1.position, n2.position, n3.position)
        ev = eigvals(Symmetric(Matrix(K)))
        @test count(<(1e-9 * maximum(ev)), ev) == 6
        # allocation contract of the in-place assembly path (Julia 1.10 heap-allocates
        # the large static matrices passed between non-inlined calls)
        X = (n1.position, n2.position, n3.position)
        shell_stiffness(sec, X...)
        @test (@allocated shell_stiffness(sec, X...)) == 0 skip = VERSION < v"1.11"
    end

    @testset "membrane: Cook's membrane converges to 23.96" begin
        v = [_cook_tri(n) for n in (4, 8, 16, 32)]
        @test all(diff(v) .> 0)                         # monotone from below
        @test v[end] ≈ 23.96 rtol = 0.01
        @test v[2] ≈ 23.96 rtol = 0.06                  # ANDES: accurate on coarse meshes
    end

    @testset "plate: square plates (thin), uniform load" begin
        @test _plate_tri(16; edge = :ss) ≈ 0.00406 rtol = 0.02
        @test _plate_tri(16; edge = :clamped) ≈ 0.00126 rtol = 0.03
    end

    @testset "shell: Scordelis–Lo roof" begin
        w = _scordelis_tri(16)
        @test w ≈ 0.3024 rtol = 0.03
    end

    @testset "patch tests on an irregular triangle mesh" begin
        # square 2 × 2, interior node moved off-grid; uniform tension along x and a
        # uniform twisting/bending state must be reproduced exactly in every element
        sec = ShellSection(Material(1000.0, 1.0, 0.25), 0.1)
        # membrane only: out-of-plane DOFs fixed; left edge holds ux (uy at one node)
        free(x, y) = x == 0 ? Bool[0, y == 0 ? 0 : 1, 0, 0, 0, 1] : Bool[1, 1, 0, 0, 0, 1]
        P = [Node([x, y, 0.0], free(x, y)) for y in (0.0, 1.0, 2.0) for x in (0.0, 1.0, 2.0)]
        P[5] = Node([0.83, 1.21, 0.0], free(0.83, 1.21))
        g(i, j) = P[3j+i+1]
        tris = AbstractElement{Float64}[]
        for j in 0:1, i in 0:1
            q = (g(i, j), g(i + 1, j), g(i + 1, j + 1), g(i, j + 1))
            append!(tris, _tri_split([ShellElement(q, sec)]; alt = isodd(i + j)))
        end
        Nx = 7.0 * 0.1                                   # σxx = 7 on t = 0.1
        # the consistent load of a constant stress state is Σₑ Lₑ·N (ANDES force
        # lumping, α_b = 3/2): boundary nodes get forces AND drilling moments
        # (∓α_b·p·ℓ²/12 per edge); interior contributions cancel
        acc = Dict{Node{Float64},Vector{Float64}}()
        for el in tris
            R, xy, _ = shell_frame(map(n -> n.position, el.nodes)...)
            Q = R[1:2, 1:2]                              # global → local in-plane
            Nl = Q * [Nx 0; 0 0] * Q'
            f = Asap._andes_parts(xy)[1] * SVector(Nl[1, 1], Nl[2, 2], Nl[1, 2])
            for (k, n) in enumerate(el.nodes)
                fg = Q' * f[3k-2:3k-1]
                acc[n] = get(acc, n, zeros(3)) .+ [fg[1], fg[2], f[3k]]
            end
        end
        loads = AbstractLoad{Float64}[]
        for (n, f) in acc
            push!(loads, NodeForce(n, [f[1], f[2], 0.0]))
            push!(loads, NodeMoment(n, [0.0, 0.0, f[3]]))
        end
        m = Model(P, tris, loads)
        solve!(m)
        for el in tris
            r = shell_resultants(m, el)
            Rl = r.frame
            # rotate the local membrane force tensor back to global: should be (Nx, 0, 0)
            Nl = [r.N[1] r.N[3]; r.N[3] r.N[2]]
            Ng = Rl[1:2, 1:2]' * Nl * Rl[1:2, 1:2]
            @test Ng ≈ [Nx 0; 0 0] atol = 1e-9
            # every point of the element, not only the centroid
            for at in ((0.0, 0.0), (1.0, 0.0), (0.2, 0.7))
                Na = shell_resultants(m, el; at = at).N
                @test [Na[1] Na[3]; Na[3] Na[2]] ≈ Nl atol = 1e-9
            end
        end
    end

    @testset "mixed quad/triangle mesh: equilibrium, area loads, resultants" begin
        sec = ShellSection(Material(3.0e4, 2.4e-3, 0.2), 0.2)
        grid, quads = grid_shell_mesh(4, 3, (i, j) -> (1.5i, 0.0, 1.0j), sec;
            fix = (i, j) -> j == 0 ? Bool[0, 0, 0, 0, 0, 0] : nothing)
        # right half triangulated
        els = AbstractElement{Float64}[]
        for el in quads
            c = sum(n -> n.position, el.nodes) / 4
            c[1] > 3.0 ? append!(els, _tri_split([el])) : push!(els, el)
        end
        loads = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81]) for el in els]...,
            [AreaLoad(el, 0.05) for el in els]...,                          # pressure along the normal
            NodeForce(grid[5, 4], [2.0, 0.5, 0.0])]
        m = Model(vec(grid), els, loads)
        solve!(m)
        # global equilibrium: reactions balance every applied load
        W = sum(el -> areal_mass(el.section) * 9.81 * Asap.area(el), els)
        Fz = -W
        Fy = sum(el -> 0.05 * Asap.area(el) * local_frame(el)[3, 2], els)
        Rsum = sum(n -> reaction(m.results, n)[1:3], grid[:, 1])
        @test Rsum ≈ -[2.0, 0.5 + Fy, Fz] rtol = 1e-9
        # section cut through the triangle zone balances the loads above it
        top = filter(el -> minimum(n -> n.position[3], el.nodes) >= 1.0 - 1e-9, els)
        cut = filter(n -> n.position[3] ≈ 1.0, vec(grid))
        F, _ = section_cut(m, top, cut)
        Wtop = sum(el -> areal_mass(el.section) * 9.81 * Asap.area(el), top)
        Ftop = sum(el -> 0.05 * Asap.area(el) * local_frame(el)[3, 2], top)
        @test F ≈ -[2.0, 0.5 + Ftop, -Wtop] rtol = 1e-8
        # ShellGeo on a mixed mesh
        g = ShellGeo(m)
        ntri = count(el -> el isa ShellElement{Float64,3}, els)
        nquad = length(els) - ntri
        @test size(shell_faces(g), 2) == ntri + 2nquad
        @test length(shell_face_owner(g)) == size(shell_faces(g), 2)
        @test length(nodal_average(g, g.N1)) == length(m.nodes)
        @test all(isfinite, g.M) && all(isfinite, g.Q)
        # pure path ≡ in-place (design-dependent triangle loads)
        st = extract_state(m)
        @test Asap.assemble_Pf(m.cache, st) ≈ m.cache.Pf rtol = 1e-12
        @test compliance(m, st) ≈ m.results.compliance rtol = 1e-10
        @test solve(m, st) ≈ m.results.u rtol = 1e-9
        # load cases and envelopes on triangles
        cr = solve_cases!(m)
        tri = first(filter(el -> el isa ShellElement{Float64,3}, els))
        env = shell_envelope(m, tri, cr, [LoadCombination(:c1, [:LC1 => 1.2]), LoadCombination(:c2, [:LC1 => 0.9])])
        @test all(env.Nmax .>= env.Nmin)
    end
end
