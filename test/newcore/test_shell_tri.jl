# S7b/S7c: triangle kernels against their primary sources.
#   ANDES-OPT membrane: Felippa (2003), CMAME 192:2125 — worked example (eqs 50–52),
#     β₀ (eqs 42, 44), exact pure bending of rectangular mesh units (Section 5).
#   MITC3+ plate: Lee, Lee & Bathe (2014), Comput. Struct. 138:12 — single-element
#     eigenvalues (Table 3, d = 1/10 000), plus patch and isotropy checks.

using LinearAlgebra, StaticArrays

_iso_plane(E, ν) = E / (1 - ν^2) * SMatrix{3,3}(1, ν, 0, ν, 1, 0, 0, 0, (1 - ν) / 2)
_centered(xy) = (c = (xy[1] + xy[2] + xy[3]) / 3; map(p -> p - c, xy))

@testset "Triangle kernels (S7b, S7c)" begin

    @testset "ANDES-OPT membrane: Felippa (2003) worked example" begin
        xy = (SVector(0.0, 0.0), SVector(102 / 25, -86 / 25), SVector(85 / 25, 285 / 250))
        A = (1 / 8) * _iso_plane(120.0, 0.25)
        K = andes_membrane_stiffness(xy, A)
        ev = sort(eigvals(Symmetric(Matrix(K))); rev = true)
        @test ev[1:6] ≈ [52.913, 43.834, 26.434, 11.181, 1.8722, 0.64900] rtol = 1e-4
        @test maximum(abs, ev[7:9]) < 1e-12 * ev[1]                 # exactly 3 rigid modes
        @test K ≈ K' rtol = 1e-14
        # selected entries of K_R (eq 52)
        @test K[1, 1] ≈ 10.392 rtol = 1e-4
        @test K[9, 9] ≈ 34.716 rtol = 1e-4
        @test K[3, 9] ≈ -0.13702 rtol = 1e-4
        @test K[5, 6] ≈ -5.5238 rtol = 1e-4
        # translation invariance: the kernel depends on coordinate differences only
        @test andes_membrane_stiffness(map(p -> p + SVector(3.0, -7.0), xy), A) ≈ K rtol = 1e-12
    end

    @testset "ANDES-OPT β₀" begin
        for ν in (0.0, 0.2, 0.25, 0.3)
            @test andes_beta0(_iso_plane(1.0, ν)) ≈ (1 - 4ν^2) / 2 rtol = 1e-12
        end
        @test andes_beta0(_iso_plane(1.0, 0.499)) == 0.01          # ½(1−4ν²) = 0.002: clamped for rank sufficiency
        # orthotropic (cracked-wall modifiers): invariant under in-plane rotation of axes
        Ao = SMatrix{3,3}(1.0, 0.1, 0, 0.1, 0.35, 0, 0, 0, 0.2)
        c, s = cos(0.4), sin(0.4)
        T = SMatrix{3,3}(c^2, s^2, -2c * s, s^2, c^2, 2c * s, c * s, -c * s, c^2 - s^2)   # strain rotation
        Ar = inv(T)' * Ao * inv(T)
        @test andes_beta0(Ar) ≈ andes_beta0(Ao) rtol = 1e-10
        @test andes_beta0(Ao) > 0
    end

    @testset "ANDES-OPT: constant-strain patch test and rigid motions" begin
        xy = _centered((SVector(0.3, -0.2), SVector(2.1, 0.4), SVector(0.8, 1.7)))
        A = 0.2 * _iso_plane(3.0e4, 0.2)
        K = andes_membrane_stiffness(xy, A)
        Ar = abs((xy[2] - xy[1])[1] * (xy[3] - xy[1])[2] - (xy[2] - xy[1])[2] * (xy[3] - xy[1])[1]) / 2
        for ε in (SVector(1e-3, 0, 0), SVector(0, 2e-3, 0), SVector(0, 0, 1.5e-3), SVector(1e-3, -4e-4, 7e-4))
            # u = εxx x + ½γ y, v = ½γ x + εyy y: zero rotation, so ψ = 0
            d = vcat(map(p -> SVector(ε[1] * p[1] + ε[3] / 2 * p[2], ε[3] / 2 * p[1] + ε[2] * p[2], 0.0), xy)...)
            @test dot(d, K * d) / 2 ≈ Ar * dot(ε, A * ε) / 2 rtol = 1e-12
        end
        # rigid rotation: u = −ωy, v = ωx, ψ = ω
        ω = 1e-3
        d = vcat(map(p -> SVector(-ω * p[2], ω * p[1], ω), xy)...)
        @test norm(K * d) < 1e-12 * norm(K) * ω
    end

    @testset "ANDES-OPT: exact in-plane pure bending of rectangular mesh units" begin
        # Felippa Section 5.1: a beam strip under constant moment, modelled by one
        # rectangle a × b split into two triangles. The energy ratio element/exact is 1
        # for every aspect ratio and Poisson ratio (the optimality property).
        for ν in (0.0, 0.25), γ in (1.0, 2.0, 4.0, 16.0)
            b = 1.0; a = γ * b; κ = 1e-3
            A = _iso_plane(1.0, ν)
            P = (SVector(0.0, -b / 2), SVector(a, -b / 2), SVector(a, b / 2), SVector(0.0, b / 2))
            # exact field (Felippa eq 35): u = −κxy, v = ½κ(x² + νy²); drilling ψ = κx
            dof(p) = SVector(-κ * p[1] * p[2], κ * (p[1]^2 + ν * p[2]^2) / 2, κ * p[1])
            U = sum(((1, 2, 3), (1, 3, 4))) do tri
                xy = map(i -> P[i], tri)
                K = andes_membrane_stiffness(xy, A)
                d = vcat(map(dof, xy)...)
                dot(d, K * d) / 2
            end
            Uexact = a * b^3 / 12 * κ^2 / 2          # σ = −κy (E = 1, plane stress), per unit thickness
            @test U / Uexact ≈ 1 rtol = 1e-10
        end
    end

    @testset "MITC3+ plate: Lee, Lee & Bathe (2014) Table 3" begin
        E, ν, t = 1.7472e7, 0.3, 1e-4
        D = t^3 / 12 * _iso_plane(E, ν)
        S = E / (2(1 + ν)) * t * SMatrix{2,2}(1.0, 0, 0, 1)   # continuum shell: no shear correction
        xy = _centered((SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(0.0, 1.0)))
        ev = sort(eigvals(Symmetric(Matrix(Asap._mitc3p_K11(xy, D, S)))))
        @test maximum(abs, ev[1:3]) < 1e-9                        # 3 rigid plate modes
        # B, B, B, B+, T (in-plane twisting), B+, S, S — 5 significant digits as published
        @test ev[4:11] ≈ [6.6685e-7, 7.9621e-7, 2.4921e-6, 8.3107e-6, 1.3599e-5, 1.4128e-5,
                          4.6667e2, 1.1760e3] rtol = 1e-4
    end

    @testset "MITC3+ plate: condensation, isotropy, constant-curvature patch test" begin
        D = (0.2^3 / 12) * _iso_plane(3.0e4, 0.2)
        S = (5 / 6) * 0.2 * (3.0e4 / 2.4) * SMatrix{2,2}(1.0, 0, 0, 1)
        xy = _centered((SVector(0.3, -0.2), SVector(2.1, 0.4), SVector(0.8, 1.7)))
        K = mitc3p_plate_stiffness(xy, D, S)
        ev = sort(eigvals(Symmetric(Matrix(K))))
        @test maximum(abs, ev[1:3]) < 1e-10 * ev[end]            # w, θx, θy rigid modes only
        @test ev[4] > 1e-8 * ev[end]
        # spatial isotropy: a cyclic relabelling of the nodes permutes K, nothing else
        xyc = (xy[2], xy[3], xy[1])
        Kc = mitc3p_plate_stiffness(xyc, D, S)
        p = SVector(4, 5, 6, 7, 8, 9, 1, 2, 3)
        @test Kc ≈ K[p, p] rtol = 1e-10
        # constant curvature with zero shear: β = −∇w, w = −½(kx x² + ky y² + kxy xy)
        for k in (SVector(1e-3, 0, 0), SVector(0, 1e-3, 0), SVector(0, 0, 1e-3), SVector(4e-4, -2e-4, 3e-4))
            w(p) = -(k[1] * p[1]^2 + k[2] * p[2]^2 + k[3] * p[1] * p[2]) / 2
            βx(p) = k[1] * p[1] + k[3] * p[2] / 2      # = −∂w/∂x
            βy(p) = k[2] * p[2] + k[3] * p[1] / 2      # = −∂w/∂y
            # (w, θx, θy) with βx = θy, βy = −θx
            d = vcat(map(p -> SVector(w(p), -βy(p), βx(p)), xy)...)
            Ar = abs((xy[2] - xy[1])[1] * (xy[3] - xy[1])[2] - (xy[2] - xy[1])[2] * (xy[3] - xy[1])[1]) / 2
            @test dot(d, K * d) / 2 ≈ Ar * dot(k, D * k) / 2 rtol = 1e-10
        end
    end
end
