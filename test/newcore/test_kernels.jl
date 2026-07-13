# Element kernel validation against the pinned Phase 0 oracles: the general
# fixity-factor stiffness must reproduce every legacy closed-form release
# matrix EXACTLY at the k ∈ {0, ∞} limits, and the transformation must match
# the legacy R!/LCS on the skewed characterization element.

# fixture geometry/section (see char_fef_element in characterization/models.jl)
const FEF_X1 = [0.0, 0.0, 0.0]
const FEF_X2 = [3000.0, 1000.0, 2000.0]
const FEF_SECTION_LEGACY = (A=1e4, E=200.0, G=77.0, Ix=8e7, Iy=3e7, J=5e6)

fef_section() = AsapNext.Section(
    AsapNext.Material(FEF_SECTION_LEGACY.E, FEF_SECTION_LEGACY.G, 1.0, 0.3),
    FEF_SECTION_LEGACY.A, FEF_SECTION_LEGACY.Ix, FEF_SECTION_LEGACY.Iy, FEF_SECTION_LEGACY.J)

# function barrier so @allocated measures the kernel, not global boxing
function _frame_kernel_allocs(sec, ends, x1, x2, Ψ)
    AsapNext.frame_stiffness(sec, ends, x1, x2, Ψ)  # warm-up/compile
    return @allocated AsapNext.frame_stiffness(sec, ends, x1, x2, Ψ)
end

@testset "Element kernels vs oracles" begin
    N = AsapNext
    sec = fef_section()
    L = N.element_length(FEF_X1, FEF_X2)
    Ψ = pi / 2   # legacy constructor default

    @testset "geometry & transformation" begin
        @test L ≈ FIXTURES["fef_fixedfixed/length"]

        Λ = N.local_frame(FEF_X1, FEF_X2, Ψ)
        R_legacy = FIXTURES["fef_fixedfixed/R"]
        @test Matrix(Λ) ≈ R_legacy[1:3, 1:3] rtol = 1e-12

        # Λ rows are the legacy LCS vectors
        lcs = FIXTURES["fef_fixedfixed/LCS"]
        for i in 1:3
            @test collect(Λ[i, :]) ≈ lcs[i] rtol = 1e-12
        end

        # right-handed orthonormal triad
        @test Λ * Λ' ≈ I(3) atol = 1e-14
        @test det(Λ) ≈ 1.0

        # vertical member (parallel to global Y) hits the special-case branch
        Λv = N.local_frame([0.0, 0.0, 0.0], [0.0, 5.0, 0.0], 0.0)
        @test Λv * Λv' ≈ I(3) atol = 1e-14
        @test collect(Λv[1, :]) ≈ [0.0, 1.0, 0.0]
    end

    @testset "local stiffness ≡ legacy closed forms: $sym" for sym in
                                                               (:fixedfixed, :fixedfree, :freefixed, :freefree, :joist)
        k = N.local_stiffness(sec, L, N.EndConditions(sym))
        @test Matrix(k) ≈ FIXTURES["fef_$sym/local_K"] rtol = 1e-12

        # global stiffness through the blockwise transform matches the pinned
        # legacy R' * k * R
        Kg = N.frame_stiffness(sec, N.EndConditions(sym), FEF_X1, FEF_X2, Ψ)
        @test Matrix(Kg) ≈ FIXTURES["fef_$sym/global_K"] rtol = 1e-10

        # symmetry of both
        @test Matrix(k) ≈ Matrix(k') rtol = 1e-12
        @test Matrix(Kg) ≈ Matrix(Kg') rtol = 1e-10
    end

    @testset "semi-rigid interpolates between limits" begin
        rigid = N.local_stiffness(sec, L, N.EndConditions(:fixedfixed))
        hinged = N.local_stiffness(sec, L, N.EndConditions(:fixedfree))

        # a very stiff end spring ≈ rigid; a very soft one ≈ hinged
        soft = N.EndConditions(N.rigid_end(), N.EndSprings(Inf, 1e-8, 1e-8, 1e-8))
        stiff = N.EndConditions(N.rigid_end(), N.EndSprings(Inf, 1e12, 1e12, 1e12))
        @test Matrix(N.local_stiffness(sec, L, soft)) ≈ Matrix(hinged) rtol = 1e-6
        @test Matrix(N.local_stiffness(sec, L, stiff)) ≈ Matrix(rigid) rtol = 1e-4

        # intermediate spring: bending diagonal strictly between the limits
        semi = N.EndConditions(N.rigid_end(), N.EndSprings(Inf, Inf, 1e9, 1e9))
        ks = N.local_stiffness(sec, L, semi)
        @test hinged[2, 2] < ks[2, 2] < rigid[2, 2]

        # fixity factor limits
        @test N.fixity_factor(Inf, 1e9, 10.0) == 1.0
        @test N.fixity_factor(0.0, 1e9, 10.0) == 0.0
        @test 0 < N.fixity_factor(1e8, 1e9, 10.0) < 1
    end

    @testset "truss kernel" begin
        n = normalize(FEF_X2 - FEF_X1)
        Kt = N.truss_stiffness(sec, FEF_X1, FEF_X2)
        ea_L = N.EA(sec) / L
        @test Matrix(Kt)[1:3, 1:3] ≈ ea_L * (n * n') rtol = 1e-12
        @test Matrix(Kt)[1:3, 4:6] ≈ -ea_L * (n * n') rtol = 1e-12
        @test Matrix(Kt) ≈ Matrix(Kt') rtol = 1e-12

        # rigid-body translation produces zero force
        u_rigid = [1.0, 2.0, 3.0, 1.0, 2.0, 3.0]
        @test norm(Matrix(Kt) * u_rigid) < 1e-8 * ea_L
    end

    @testset "kernels are allocation-free and generic" begin
        ends = N.EndConditions(:fixedfixed)
        x1s, x2s = SVector{3}(FEF_X1), SVector{3}(FEF_X2)
        @test _frame_kernel_allocs(sec, ends, x1s, x2s, Ψ) == 0

        # scalar-generic: BigFloat flows through the whole kernel
        mb = N.Material(big"200.0", big"77.0", big"1.0", big"0.3")
        sb = N.Section(mb, big"1e4", big"8e7", big"3e7", big"5e6")
        kb = N.local_stiffness(sb, big(L), N.EndConditions(:fixedfixed))
        @test eltype(kb) == BigFloat
        @test Float64.(Matrix(kb)) ≈ FIXTURES["fef_fixedfixed/local_K"] rtol = 1e-12
    end
end
