# Shell S4: AreaLoad / SelfWeight on shells, stress resultants, section cuts.
# Oracles: statics (exact for consistent loads and nodal-force cuts), the
# uniform-tension membrane state, Timoshenko & Woinowsky-Krieger plate moment.

using StaticArrays, LinearAlgebra
include("shell_test_utils.jl")

# area centroid of a planar quad (two triangles)
function _quad_centroid(X)
    a1 = norm(cross(X[2] - X[1], X[3] - X[1])) / 2
    a2 = norm(cross(X[3] - X[1], X[4] - X[1])) / 2
    return (a1 * (X[1] + X[2] + X[3]) / 3 + a2 * (X[1] + X[3] + X[4]) / 3) / (a1 + a2)
end

# vertical cantilever wall in the XZ plane: length B (x), height H (z), clamped base
function _wall(nx, nz; B = 4.0, H = 6.0, t = 0.2, E = 3.0e7, ν = 0.2, ρ = 2.4)
    sec = ShellSection(Material(E, ρ, ν), t)
    grid, els = shell_mesh(nx, nz, (i, j) -> (B * i / nx, 0.0, H * j / nz), sec;
        fix = (i, j) -> j == 0 ? falses(6) |> collect : nothing)
    return grid, els, sec
end

@testset "Shell loads, resultants, section cuts (S4)" begin

    @testset "consistent area loads: totals and moment balance" begin
        sec = ShellSection(Material(100.0, 1.0, 0.3), 0.1)
        X = (SVector(0.0, 0, 0), SVector(2.0, 0.3, 0.1), SVector(2.2, 1.9, 0.4), SVector(-0.1, 1.5, 0.3))
        X = (X[1], X[2], X[3], X[1] + 0.9 * (X[3] - X[2]) + 0.05 * (X[2] - X[1]))   # flat, irregular quad
        n = [Node(collect(x), :free) for x in X]
        el = ShellElement(n, sec)
        q = SVector(0.3, -0.7, 1.1)
        q̃, R = Asap.shell_load_vector(sec, X..., q)
        fg = Asap._rotate24(R', q̃)                                   # global FEF = −consistent load
        F = -sum(SVector(fg[6k-5], fg[6k-4], fg[6k-3]) for k in 1:4)
        @test F ≈ q * Asap.area(el) rtol = 1e-12
        c = sum(X) / 4
        Mc = -sum(cross(X[k] - c, SVector(fg[6k-5], fg[6k-4], fg[6k-3])) +
                  SVector(fg[6k-2], fg[6k-1], fg[6k]) for k in 1:4)
        # the resultant of a uniform traction acts at the area centroid
        @test norm(Mc - cross(_quad_centroid(X) - c, q * Asap.area(el))) < 1e-10 * norm(q) * Asap.area(el)
        @test_throws ArgumentError AreaLoad(FrameElement(n[1], n[2], Section(Material(1.0, 1.0, 0.3), 1.0, 1.0, 1.0, 1.0)), [0.0, 0, 1])
    end

    @testset "SelfWeight on a wall: base reaction and storey cut are exact" begin
        grid, els, sec = _wall(4, 6)
        g = SVector(0.0, 0.0, -9.81)
        m = Model(vec(grid), els, AbstractLoad{Float64}[SelfWeight(el; g = g) for el in els])
        solve!(m)
        W = areal_mass(sec) * 9.81 * 4.0 * 6.0
        Rz = sum(reaction(m.results, grid[i, 1])[3] for i in 1:5)
        @test Rz ≈ W rtol = 1e-10
        # cut at storey level z = 3: the part below carries the weight above
        below = [el for el in els if all(n -> n.position[3] <= 3.0 + 1e-9, el.nodes)]
        cut = [grid[i, 4] for i in 1:5]
        F, M = section_cut(m, below, cut)
        @test F[3] ≈ -W / 2 rtol = 1e-10          # the upper half pushes down on the lower half
        @test abs(F[1]) < 1e-9 * W
        @test abs(M[2]) < 1e-9 * W                 # symmetric: no overturning
    end

    @testset "section cut: lateral load on a cantilever wall (exact statics)" begin
        grid, els, sec = _wall(4, 6)
        P = 100.0
        m = Model(vec(grid), els, AbstractLoad{Float64}[NodeForce(grid[i, 7], [P / 5, 0.0, 0.0]) for i in 1:5])
        solve!(m)
        for (jc, zc) in ((2, 1.0), (4, 3.0), (6, 5.0))
            below = [el for el in els if all(n -> n.position[3] <= zc + 1e-9, el.nodes)]
            cut = [grid[i, jc] for i in 1:5]
            F, M = section_cut(m, below, cut; about = [2.0, 0.0, zc])
            @test F[1] ≈ P rtol = 1e-10            # storey shear
            @test abs(F[3]) < 1e-9 * P             # no axial
            @test M[2] ≈ P * (6.0 - zc) rtol = 1e-10   # overturning about y
            # pier axes: rows (along wall, through thickness, vertical)
            Fp, Mp = section_cut(m, below, cut; about = [2.0, 0.0, zc], axes = [1.0 0 0; 0 1 0; 0 0 1])
            @test Fp ≈ F
        end
        # base reactions balance the load
        @test sum(reaction(m.results, grid[i, 1])[1] for i in 1:5) ≈ -P rtol = 1e-10
    end

    @testset "membrane resultants: uniform tension is reproduced exactly" begin
        sec = ShellSection(Material(1000.0, 1.0, 0.3), 0.5)
        # distorted strip, uniform edge tension σ·t = 2 [force/length] on x = 4
        xy(i, j) = (i * 1.0 + (0 < i < 4 && j == 1 ? 0.2 : 0.0), 2.0 * j / 2 + (i == 2 && j == 1 ? 0.15 : 0.0))
        # ψ restrained on the supported edge: the exact field has ψ = 0, and a
        # free ψ there would let the Allman edge term bulge the "support"
        grid, els = membrane_mesh(4, 2, xy, sec; fix = (i, j) -> i == 0 ? (true, j == 0, true) : nothing)
        Nt = 2.0
        # consistent edge load of a uniform normal traction on the Allman element:
        # forces Nt·l/2 per edge end PLUS drilling moments ∓Nt·l²/12 (edge from
        # J to K: −at J, +at K); interior node moments cancel
        loads = AbstractLoad{Float64}[[NodeForce(grid[5, j+1], [Nt * (j == 1 ? 1.0 : 0.5), 0.0, 0.0]) for j in 0:2]...,
                 NodeMoment(grid[5, 1], [0.0, 0.0, -Nt / 12]), NodeMoment(grid[5, 3], [0.0, 0.0, Nt / 12])]
        m = Model(vec(grid), els, AbstractLoad{Float64}[loads...])
        solve!(m)
        for el in els, at in ((0.0, 0.0), (0.5, -0.7))
            r = shell_resultants(m, el; at = at)
            # local axes = global for these XY elements up to the edge-1 direction: rotate back
            Rp = r.frame[1:2, 1:2]
            Nmat = Rp' * [r.N[1] r.N[3]; r.N[3] r.N[2]] * Rp
            @test Nmat ≈ [Nt 0; 0 0] atol = 1e-9
            @test norm(r.M) < 1e-12 && norm(r.Q) < 1e-12
        end
    end

    @testset "plate resultants: SS square plate centre moment" begin
        # Timoshenko & Woinowsky-Krieger: Mx = My = 0.0479 q a² at the centre (ν = 0.3)
        n = 16
        sec = ShellSection(Material(1.0e4, 1.0, 0.3), 0.01)
        onedge(i, j) = i == 0 || j == 0 || i == n || j == n
        grid, els = shell_mesh(n, n, (i, j) -> (i / n, j / n, 0.0), sec;
            fix = (i, j) -> onedge(i, j) ? Bool[0, 0, 0, 1, 1, 0] : Bool[0, 0, 1, 1, 1, 0])
        m = Model(vec(grid), els, AbstractLoad{Float64}[AreaLoad(el, [0.0, 0.0, -1.0]) for el in els])
        solve!(m)
        # AreaLoad ≡ lumped loads on rectangles: same deflection as the S3 helper
        D = 1.0e4 * 0.01^3 / (12 * (1 - 0.09))
        @test -displacement(m.results, grid[n÷2+1, n÷2+1])[3] * D ≈ square_plate(n) rtol = 1e-10
        # element whose top-right corner is the plate centre
        el = only(filter(e -> e.nodes[3] === grid[n÷2+1, n÷2+1], els))
        r = shell_resultants(m, el; at = (1.0, 1.0))
        @test abs(r.M[1]) ≈ 0.0479 rtol = 0.03
        @test abs(r.M[2]) ≈ 0.0479 rtol = 0.03
        @test r.M[1] < 0                            # sagging: +e₃ (top) face in compression
        # equilibrium of a quarter: edge shears carry the load (centre element Q ≈ 0 by symmetry)
        @test norm(r.Q) < 0.05
    end

    @testset "cases and combinations" begin
        grid, els, sec = _wall(2, 3)
        loads = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81], case = :D) for el in els]...,
                                      [AreaLoad(el, [0.0, 1.5, 0.0]; case = :W) for el in els]...]
        m = Model(vec(grid), els, loads)
        cr = solve_cases!(m)
        combo = LoadCombination(:ULS, [:D => 1.2, :W => 1.6])
        rc = combine(cr, combo)
        # the same loads, factored, on an independent model
        g2, e2, _ = _wall(2, 3)
        l2 = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81], factor = 1.2) for el in e2]...,
                                   [AreaLoad(el, [0.0, 1.6 * 1.5, 0.0]) for el in e2]...]
        m2 = Model(vec(g2), e2, l2)
        solve!(m2)
        @test rc.u ≈ m2.results.u rtol = 1e-10
        for k in eachindex(els)
            ra = shell_resultants(m, els[k]; results = rc, factors = Dict(:D => 1.2, :W => 1.6))
            rb = shell_resultants(m2, e2[k])
            @test ra.N ≈ rb.N rtol = 1e-9 atol = 1e-9
            @test ra.M ≈ rb.M rtol = 1e-9 atol = 1e-9
        end
    end
end
