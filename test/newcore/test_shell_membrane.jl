# Shell S2: ShellSection + the drilling membrane — Ibrahimbegović–Taylor–Wilson
# (ITW 1990, D-type) with two complementary Taylor-corrected incompatible modes
# (default) or as published (`enhanced = false`).
#
# Oracles: exact solutions (pure bending, patch test) and fine-mesh references
# computed through this pipeline AND cross-checked against an independent
# conforming Q4 code (slender cantilever: 256×64 Q4 → 0.35601).

using StaticArrays, LinearAlgebra
include("shell_test_utils.jl")

const CANT_REF = 0.3559997546941379      # 128×32, l/h = 4 (Q4 256×64: 0.35601)
const SQUAT_REF = 0.009155330838561265   # 64×64, l/h = 1
const COOK_REF = 23.96                   # right-edge midpoint (literature)

@testset "Shell membrane (S2)" begin

    @testset "ShellSection rigidities" begin
        mat = Material(30.0, 1.0, 0.25)
        s = ShellSection(mat, 0.5)
        c = 30.0 / (1 - 0.25^2)
        @test membrane_rigidity(s) ≈ 0.5 * [c 0.25c 0; 0.25c c 0; 0 0 mat.G]
        @test bending_rigidity(s) ≈ (0.5^3 / 12) * [c 0.25c 0; 0.25c c 0; 0 0 mat.G]
        @test shear_rigidity(s) ≈ (5 / 6) * mat.G * 0.5 * I(2)
        @test drilling_rigidity(s) ≈ mat.G * 0.5
        @test areal_mass(s) ≈ 0.5
        @test thickness(s) == 0.5
        @test_throws ArgumentError ShellSection(mat, 0.0)
        m = ShellRigiditySection(s; f11 = 0.5, m11 = 0.35, m22 = 0.35, m12 = 0.35, v13 = 0.8)
        @test membrane_rigidity(m)[1, 1] ≈ 0.5 * membrane_rigidity(s)[1, 1]
        @test membrane_rigidity(m)[1, 2] ≈ sqrt(0.5) * membrane_rigidity(s)[1, 2]
        @test bending_rigidity(m) ≈ 0.35 * bending_rigidity(s)
        @test shear_rigidity(m)[1, 1] ≈ 0.8 * shear_rigidity(s)[1, 1]
        @test drilling_rigidity(m) ≈ membrane_rigidity(m)[3, 3]
        @test_throws ArgumentError ShellRigiditySection(-membrane_rigidity(s),
            bending_rigidity(s), shear_rigidity(s), 1.0)
    end

    @testset "kernel: rank, rigid modes, invariance, allocations" begin
        sec = ShellSection(Material(100.0, 1.0, 0.3), 1.0)
        X = (SVector(0.0, 0, 0), SVector(2.0, 0.3, 0), SVector(2.2, 1.9, 0), SVector(-0.1, 1.5, 0))
        R, xy, warp = shell_frame(X...)
        @test warp ≈ 0 atol = 1e-14
        @test R * R' ≈ I(3)
        for enhanced in (true, false)
            Km = itw_membrane_stiffness(xy, membrane_rigidity(sec), drilling_rigidity(sec); enhanced)
            @test Km ≈ Km'
            ev = eigvals(Symmetric(Matrix(Km)))
            @test count(<(1e-10 * ev[end]), ev) == 3          # u, v, in-plane rigid rotation only
            θ = 1e-2                                           # rigid rotation incl. ψ = θ: zero energy
            a = vcat([SVector(-θ * p[2], θ * p[1], θ) for p in xy]...)
            @test norm(Km * a) < 1e-12 * norm(Km) * norm(a)
            # constant ψ with no displacement is NOT a zero-energy mode (ITW §4.1)
            ψ = vcat([SVector(0.0, 0.0, 1.0) for _ in 1:4]...)
            @test ψ' * Km * ψ > 1e-3 * ev[end]
        end
        # global kernel: rotation invariance of the spectrum
        K = shell_stiffness(sec, X...)
        Q = SMatrix{3,3}(Matrix(qr([0.3 -1.2 0.4; 0.9 0.1 -0.7; 0.2 0.5 1.1]).Q))
        Kr = shell_stiffness(sec, map(x -> Q * x + SVector(1.0, -2.0, 3.0), X)...)
        @test sort(eigvals(Symmetric(Matrix(Kr)))) ≈ sort(eigvals(Symmetric(Matrix(K)))) rtol = 1e-10 atol = 1e-10
        # allocation contract of the in-place assembly path
        shell_stiffness(sec, X...)
        # Julia 1.10 heap-allocates the large static matrices the shell kernels pass
        # between non-inlined calls (12×12 and up); 1.11+ keeps them on the stack
        @test (@allocated shell_stiffness(sec, X...)) == 0 skip = VERSION < v"1.11"
    end

    @testset "patch test: constant strain + rotation, distorted patch" begin
        # MacNeal–Harder 5-element patch, Dirichlet data from the linear field
        # u = a(x + y/2), v = a(y + x) ⇒ constant strain and rotation ω = a/4;
        # interior nodes must reproduce the field — including ψ = ω — exactly
        E, ν, a = 1e6, 0.25, 1e-3
        sec = ShellSection(Material(E, 1.0, ν), 0.001)
        P = [SVector(0.0, 0.0), SVector(0.24, 0.0), SVector(0.24, 0.12), SVector(0.0, 0.12),
             SVector(0.04, 0.02), SVector(0.18, 0.03), SVector(0.16, 0.08), SVector(0.08, 0.08)]
        quads = [(1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (4, 1, 5, 8), (5, 6, 7, 8)]
        field(p) = [a * (p[1] + p[2] / 2), a * (p[2] + p[1]), a / 4]
        U = vcat(field.(P)...)
        for enhanced in (true, false)
            Kg = zeros(24, 24)
            for q in quads
                Ke = itw_membrane_stiffness(Tuple(P[i] for i in q), membrane_rigidity(sec),
                    drilling_rigidity(sec); enhanced)
                d = vcat((3 * (i - 1) .+ (1:3) for i in q)...)
                Kg[d, d] .+= Ke
            end
            Ui = -Kg[13:24, 13:24] \ (Kg[13:24, 1:12] * U[1:12])
            @test Ui ≈ U[13:24] rtol = 1e-10
        end
    end

    @testset "pure bending, regular mesh (ITW Table 1 LC1)" begin
        # E = 100, ν = 0.3, l = 10, h = 1, one row of 6 elements; pin at the
        # bottom-left, roller at the bottom-right, end couples P·h = 1 from
        # nodal forces; beam theory = plane-stress elasticity here: midspan
        # deflection ML²/8EI = 1.5, end rotation ML/2EI = 0.6
        sec = ShellSection(Material(100.0, 1.0, 0.3), 1.0)
        fix(i, j) = (i, j) == (0, 0) ? (true, true, false) : (i, j) == (6, 0) ? (false, true, false) : nothing
        grid, els = membrane_mesh(6, 1, (i, j) -> (10i / 6, Float64(j)), sec; fix = fix)
        loads = [NodeForce(grid[1, 2], [-1.0, 0, 0]), NodeForce(grid[1, 1], [1.0, 0, 0]),
                 NodeForce(grid[7, 2], [1.0, 0, 0]), NodeForce(grid[7, 1], [-1.0, 0, 0])]
        m = Model(vec(grid), els, AbstractLoad{Float64}[loads...])
        solve!(m)
        @test displacement(m.results, grid[4, 1])[2] ≈ 1.5 rtol = 1e-10       # exact (enhanced)
        @test displacement(m.results, grid[4, 2])[2] ≈ 1.5 rtol = 1e-10
        @test abs(displacement(m.results, grid[7, 1])[6]) ≈ 0.6 rtol = 1e-10  # drilling DOF = rotation
        # self-equilibrated load ⇒ zero vertical support reactions (Asap's
        # convention: a load applied AT a support is not subtracted)
        @test abs(reaction(m.results, grid[1, 1])[2]) < 1e-10
        @test abs(reaction(m.results, grid[7, 1])[2]) < 1e-10
    end

    @testset "cantilever: coarse-mesh accuracy vs converged reference" begin
        @test cantilever_wall(4, 1) / CANT_REF > 0.975
        @test cantilever_wall(8, 2) / CANT_REF ≈ 1 atol = 0.005
        @test cantilever_wall(16, 4) / CANT_REF ≈ 1 atol = 0.005
        # squat (l/h = 1) wall — the shear-wall case — two elements per side
        @test cantilever_wall(2, 2; L = 12.0, H = 12.0) / SQUAT_REF ≈ 1 atol = 0.04
        @test cantilever_wall(4, 4; L = 12.0, H = 12.0) / SQUAT_REF ≈ 1 atol = 0.02
    end

    @testset "Cook's membrane (midpoint reference 23.96)" begin
        c2, c4, c8 = cook_mid(2), cook_mid(4), cook_mid(8)
        @test c2 < c4 < c8 < COOK_REF
        @test c8 / COOK_REF > 0.985
    end

    @testset "element in a model: validation and geometry" begin
        sec = ShellSection(Material(100.0, 1.0, 0.3), 0.2)
        n = [Node([0.0, 0, 0], :fixed), Node([2.0, 0, 0], :fixed), Node([2.0, 0, 1.5], :free), Node([0.0, 0, 1.5], :free)]
        el = ShellElement(n, sec, :wall)
        @test Asap.area(el) ≈ 3.0
        @test Asap.local_frame(el)[3, :] ≈ [0, -1, 0]        # XZ wall, CCW about −y
        @test_throws ArgumentError ShellElement((n[1], n[2]), sec)          # 2 nodes: not a shell
        @test_throws ArgumentError ShellElement((n[1], n[1], n[3], n[4]), sec)
        bad = ShellElement((n[1], n[3], n[2], n[4]), sec)     # bow-tie ordering
        @test_throws ArgumentError process!(Model(n, AbstractElement{Float64}[bad],
            AbstractLoad{Float64}[NodeForce(n[3], [1.0, 0, 0])]))
        m = Model(n, AbstractElement{Float64}[el], AbstractLoad{Float64}[NodeForce(n[3], [1.0, 0, 0])])
        @test volume(m) ≈ 0.2 * 3.0
        for x in (sec, ShellRigiditySection(sec; m11 = 0.35), el)
            @test !isempty(sprint(show, MIME"text/plain"(), x))
            @test !isempty(sprint(show, x))
        end
    end

    @testset "pure path parity (shell group)" begin
        sec = ShellSection(Material(30000.0, 1.0, 0.25), 1.0)
        grid, els = membrane_mesh(4, 2, (i, j) -> (12.0i, 6.0j), sec;
            fix = (i, j) -> i == 0 ? (true, true, true) : nothing)
        m = Model(vec(grid), els, AbstractLoad{Float64}[NodeForce(grid[5, 3], [0.0, 10.0, 0.0])])
        solve!(m)
        @test compliance(m, extract_state(m)) ≈ m.results.compliance rtol = 1e-12
    end
end
