# New-core Materials/Sections tests: the accessor contract and its
# equivalences. Run via test/newcore/runtests.jl.

@testset "Materials & Sections" begin
    N = AsapNext

    @testset "Material construction & promotion" begin
        m = N.Material(200e6, 77e6, 8.0, 0.3)
        @test m isa N.Material{Float64}
        @test (m.E, m.G, m.ρ, m.ν) == (200e6, 77e6, 8.0, 0.3)

        # mixed input types promote
        m2 = N.Material(200, 77, 8, 1 // 2)
        @test m2 isa N.Material{Rational{Int}}

        # isotropic G derivation: G = E / (2(1 + ν))
        m3 = N.Material(200e6, 8.0, 0.3)
        @test m3.G ≈ 200e6 / 2.6

        # parametric scalar types flow through (the point of the rewrite)
        mb = N.Material(big"200.0e6", big"8.0", big"0.3")
        @test mb isa N.Material{BigFloat}

        # legacy presets preserved
        @test N.Steel_Nmm.E == 200e3 && N.Steel_kNm.E == 200e6
    end

    @testset "Section accessors ≡ products" begin
        m = N.Material(200e6, 77e6, 8.0, 0.3)
        s = N.Section(m, 1e-2, 8e-5, 3e-5, 5e-7)

        @test N.EA(s) ≈ 200e6 * 1e-2
        @test N.EIx(s) ≈ 200e6 * 8e-5
        @test N.EIy(s) ≈ 200e6 * 3e-5
        @test N.GJ(s) ≈ 77e6 * 5e-7
        @test N.ρA(s) ≈ 8.0 * 1e-2

        # axial-only (truss) convenience: flexural/torsional rigidities zero
        bar = N.Section(m, 5e-3)
        @test N.EA(bar) ≈ 200e6 * 5e-3
        @test N.EIx(bar) == N.EIy(bar) == N.GJ(bar) == 0.0
    end

    @testset "RigiditySection ≡ equivalent Section" begin
        m = N.Material(200e6, 77e6, 8.0, 0.3)
        s = N.Section(m, 1e-2, 8e-5, 3e-5, 5e-7)
        r = N.RigiditySection(N.EA(s), N.EIx(s), N.EIy(s), N.GJ(s), N.ρA(s))

        # a RigiditySection built from a Section's products is indistinguishable
        # through the accessor contract — the whole point of the abstraction
        for f in (N.EA, N.EIx, N.EIy, N.GJ, N.ρA)
            @test f(r) == f(s)
        end

        @test r isa N.AbstractSection{Float64}
        @test s isa N.AbstractSection{Float64}
    end

    @testset "show methods render" begin
        m = N.Material(200e6, 77e6, 8.0, 0.3)
        s = N.Section(m, 1e-2, 8e-5, 3e-5, 5e-7)
        r = N.RigiditySection(1.0, 2.0, 3.0, 4.0, 5.0)

        for x in (m, s, r)
            plain = sprint(show, MIME"text/plain"(), x)
            @test !isempty(plain) && occursin("force", plain)
            @test !isempty(sprint(show, x))
        end
        # engineer-language explanations present
        @test occursin("Poisson", sprint(show, MIME"text/plain"(), m))
        @test occursin("rigidities", sprint(show, MIME"text/plain"(), s))
    end
end
