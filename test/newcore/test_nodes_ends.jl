# New-core Node and EndConditions tests, including exact compatibility of the
# fixity map with the legacy fixDict and of the release map with the legacy
# release semantics.

@testset "Node" begin
    N = AsapNext

    @testset "construction" begin
        n = N.Node([0.0, 0.0, 3.0], :pinned)
        @test n isa N.Node{Float64}
        @test n.position === SVector(0.0, 0.0, 3.0)
        @test collect(n.fixity) == [false, false, false, true, true, true]
        @test n.id == :node && n.index == 0

        n2 = N.Node([1, 2, 3], [true, true, false, true, true, true], :roller)
        @test n2 isa N.Node{Float64}      # integer positions float-promote
        @test n2.id == :roller

        @test_throws AssertionError N.Node([0.0, 0.0], :free)
        @test_throws AssertionError N.Node([0.0, 0.0, 0.0], :notafixity)
    end

    # legacy fixDict values, frozen here as the compatibility contract
    LEGACY_FIXDICT = Dict(
        :fixed => [false, false, false, false, false, false],
        :free => [true, true, true, true, true, true],
        :xfixed => [false, true, true, true, true, true],
        :yfixed => [true, false, true, true, true, true],
        :zfixed => [true, true, false, true, true, true],
        :xfree => [true, false, false, false, false, false],
        :yfree => [false, true, false, false, false, false],
        :zfree => [false, false, true, false, false, false],
        :pinned => [false, false, false, true, true, true])

    @testset "FIXITIES ≡ legacy fixDict" begin
        for (sym, legacy) in LEGACY_FIXDICT
            @test collect(N.FIXITIES[sym]) == legacy
        end
        @test length(N.FIXITIES) == length(LEGACY_FIXDICT)
    end

    @testset "fixnode! / planarize!" begin
        n = N.Node([0.0, 0.0, 0.0], :free)
        N.fixnode!(n, :fixed)
        @test all(.!n.fixity)

        n = N.Node([0.0, 0.0, 0.0], :free)
        N.planarize!(n)                     # default :XY — fixes Tz, Rx, Ry
        @test collect(n.fixity) == [true, true, false, false, false, true]

        # matches legacy planeDict conventions (values frozen here)
        for (plane, idx) in Dict(:XY => [3, 4, 5], :YZ => [1, 5, 6], :ZX => [2, 4, 6])
            n = N.Node([0.0, 0.0, 0.0], :free)
            N.planarize!(n, plane)
            expected = trues(6)
            expected[idx] .= false
            @test collect(n.fixity) == collect(expected)
        end

        ns = [N.Node([0.0, 0.0, 0.0], :free) for _ in 1:3]
        N.planarize!(ns)
        @test all(n -> !n.fixity[3], ns)
    end
end

@testset "EndSprings / EndConditions" begin
    N = AsapNext

    @testset "constructors & limits" begin
        r = N.rigid_end()
        @test all(isinf, (r.kx, r.kt, r.ky, r.kz))

        p = N.pinned_end()
        @test isinf(p.kx) && iszero(p.kt) && iszero(p.ky) && iszero(p.kz)

        @test_throws AssertionError N.EndSprings(-1.0, 0.0, 0.0, 0.0)

        semi = N.EndSprings(Inf, Inf, 5e4, 5e4)
        @test semi.ky == 5e4

        ec = N.EndConditions(semi, N.rigid_end())
        @test ec isa N.EndConditions{Float64}
    end

    @testset "release symbols round-trip to exact limits" begin
        for sym in (:fixedfixed, :fixedfree, :freefixed, :freefree, :joist)
            ec = N.EndConditions(sym)
            @test N.release_symbol(ec) == sym
        end

        # classical semantics spot checks
        ff = N.EndConditions(:fixedfree)
        @test all(isinf, (ff.e1.kx, ff.e1.kt, ff.e1.ky, ff.e1.kz))   # start rigid
        @test isinf(ff.e2.kx) && iszero(ff.e2.kt)                     # end hinged

        joist = N.EndConditions(:joist)
        @test isinf(joist.e1.kt) && iszero(joist.e1.ky)               # torsion kept, bending released

        # a semi-rigid connection is NOT a classical release
        semi = N.EndConditions(N.EndSprings(Inf, Inf, 5e4, 5e4), N.rigid_end())
        @test N.release_symbol(semi) === nothing

        @test_throws AssertionError N.EndConditions(:notarelease)
    end

    @testset "show methods render" begin
        n = N.Node([0.0, 0.0, 3.0], :pinned, :base)
        @test occursin("Tx Ty Tz", sprint(show, MIME"text/plain"(), n))

        ec = N.EndConditions(:fixedfree)
        plain = sprint(show, MIME"text/plain"(), ec)
        @test occursin("fixedfree", plain) && occursin("rigid", plain)

        semi = N.EndConditions(N.EndSprings(Inf, Inf, 5e4, 5e4), N.rigid_end())
        @test occursin("semi-rigid", sprint(show, MIME"text/plain"(), semi))
        @test occursin("torsion", sprint(show, MIME"text/plain"(), semi.e1))
    end
end
