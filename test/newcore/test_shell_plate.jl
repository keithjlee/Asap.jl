# Shell S3: MITC4 plate bending → full flat-shell element.
# Oracles: Timoshenko & Woinowsky-Krieger (plate coefficients, ν = 0.3),
# MacNeal & Harder (1985) Scordelis–Lo roof, exact constant-curvature patch.

using StaticArrays, LinearAlgebra
include("shell_test_utils.jl")

@testset "Shell plate / full shell (S3)" begin

    @testset "kernel: rank and rigid modes" begin
        sec = ShellSection(Material(100.0, 1.0, 0.3), 0.1)
        X = (SVector(0.0, 0, 0), SVector(2.0, 0.3, 0), SVector(2.2, 1.9, 0), SVector(-0.1, 1.5, 0))
        _, xy, _ = shell_frame(X...)
        Kp = Asap.mitc4_plate_stiffness(xy, bending_rigidity(sec), shear_rigidity(sec))
        @test Kp ≈ Kp'
        ev = eigvals(Symmetric(Matrix(Kp)))
        @test count(<(1e-10 * ev[end]), ev) == 3
        for (w0, tx, ty) in ((1.0, 0, 0), (0, 1.0, 0), (0, 0, 1.0))   # w = w0 + θx·y − θy·x
            v = vcat([SVector(w0 + tx * p[2] - ty * p[1], tx, ty) for p in xy]...)
            @test norm(Kp * v) < 1e-12 * norm(Kp)
        end
        # full shell: exactly 6 rigid-body modes in an arbitrary orientation
        Q = SMatrix{3,3}(Matrix(qr([0.3 -1.2 0.4; 0.9 0.1 -0.7; 0.2 0.5 1.1]).Q))
        K = shell_stiffness(sec, map(x -> Q * x + SVector(1.0, 2.0, 3.0), X)...)
        evK = eigvals(Symmetric(Matrix(K)))
        @test count(<(1e-10 * evK[end]), evK) == 6
        shell_stiffness(sec, X...)
        @test (@allocated shell_stiffness(sec, X...)) == 0
    end

    @testset "constant-curvature patch test (distorted, thin and thick)" begin
        # Kirchhoff field w = (a x² + b y² + c xy)/2, β = −∇w ⇒ θy = −w,x, θx = w,y:
        # constant curvature, zero transverse shear; interior nodes must be exact
        P = [SVector(0.0, 0.0), SVector(0.24, 0.0), SVector(0.24, 0.12), SVector(0.0, 0.12),
             SVector(0.04, 0.02), SVector(0.18, 0.03), SVector(0.16, 0.08), SVector(0.08, 0.08)]
        quads = [(1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (4, 1, 5, 8), (5, 6, 7, 8)]
        a, b, c = 1e-3, -2e-3, 1.5e-3
        field(p) = [(a * p[1]^2 + b * p[2]^2 + c * p[1] * p[2]) / 2,
                    b * p[2] + c * p[1] / 2,           # θx = w,y
                    -(a * p[1] + c * p[2] / 2)]        # θy = −w,x
        U = vcat(field.(P)...)
        for t in (0.001, 0.05)
            sec = ShellSection(Material(1e6, 1.0, 0.25), t)
            Kg = zeros(24, 24)
            for q in quads
                Ke = Asap.mitc4_plate_stiffness(Tuple(P[i] for i in q), bending_rigidity(sec), shear_rigidity(sec))
                d = vcat((3 * (i - 1) .+ (1:3) for i in q)...)
                Kg[d, d] .+= Ke
            end
            Ui = -Kg[13:24, 13:24] \ (Kg[13:24, 1:12] * U[1:12])
            @test Ui ≈ U[13:24] rtol = 1e-8
        end
    end

    @testset "square plates vs Timoshenko & Woinowsky-Krieger (ν = 0.3)" begin
        # uniform load: SS 0.00406, clamped 0.00126 (×qa⁴/D); SS centre load 0.0116 (×Pa²/D)
        @test square_plate(16) ≈ 0.00406 rtol = 0.01
        @test square_plate(16; edge = :clamped) ≈ 0.00126 rtol = 0.02
        @test square_plate(16; point = true) ≈ 0.0116 rtol = 0.02
        # monotone convergence of the SS uniform case
        s4, s8, s16 = square_plate(4), square_plate(8), square_plate(16)
        @test abs(s16 - 0.00406) < abs(s8 - 0.00406) < abs(s4 - 0.00406)
    end

    @testset "no shear locking: thin-to-thick sweep" begin
        thin = [square_plate(8; t = 1 / r) for r in (1e2, 1e3, 1e4)]
        @test thin[3] ≈ thin[2] rtol = 1e-3          # Kirchhoff limit reached, no lock-up
        @test thin[2] ≈ 0.00406 rtol = 0.02
        thick = square_plate(8; t = 0.1)             # a/t = 10: shear deformation adds flexibility
        @test thick > thin[2]
    end

    @testset "Scordelis–Lo roof (MacNeal–Harder, 0.3024)" begin
        s4, s8, s16 = scordelis_lo(4), scordelis_lo(8), scordelis_lo(16)
        @test s16 ≈ 0.3024 rtol = 0.02
        @test abs(s16 - 0.3024) < abs(s4 - 0.3024)
        @test s8 ≈ 0.3024 rtol = 0.05
    end
end
