# Element struct + interface tests: construction, dof signatures (the
# activity contract), stiffness consistency with kernels, springs.

@testset "Elements & interface" begin
    N = AsapNext
    sec = fef_section()
    n1 = N.Node(FEF_X1, :fixed)
    n2 = N.Node(FEF_X2, :free)

    @testset "FrameElement construction" begin
        el = N.FrameElement(n1, n2, sec, :girder)
        @test el isa N.FrameElement{Float64,typeof(sec)}
        @test N.nodes(el) === (n1, n2)
        @test N.ndofs(el) == 12
        @test N.n_internal_dofs(el) == 0
        @test el.rollangle ≈ pi / 2                    # legacy default
        @test N.release_symbol(el.ends) == :fixedfixed
        @test Base.length(el) ≈ FIXTURES["fef_fixedfixed/length"]

        hinged = N.FrameElement(n1, n2, sec; release=:fixedfree)
        @test N.release_symbol(hinged.ends) == :fixedfree

        semi = N.FrameElement(n1, n2, sec,
            N.EndConditions(N.EndSprings(Inf, Inf, 5e4, 5e4), N.rigid_end()), :semi)
        @test N.release_symbol(semi.ends) === nothing

        # section is swappable (design iteration) without rebuilding
        el.section = N.Section(N.Material(200.0, 77.0, 1.0, 0.3), 2e4)
        @test N.EA(el.section) ≈ 200.0 * 2e4
    end

    @testset "dof_signature drives activity" begin
        rigid = N.FrameElement(n1, n2, sec)
        @test all(N.dof_signature(rigid))

        # hinge at the far end: the far node's ENTIRE rotation block decouples
        # (all three local rotational rows are zero there); the start block
        # stays fully active because Λ mixes its kept bending rotations with
        # torsion — activity is blockwise in global coordinates
        ff = N.FrameElement(n1, n2, sec; release=:fixedfree)
        sig = N.dof_signature(ff)
        @test sig[1:3] == (true, true, true) && sig[7:9] == (true, true, true)
        @test sig[10] == sig[11] == sig[12] == false      # far rotation block off
        @test sig[4] == sig[5] == sig[6] == true          # start block on (blockwise)

        # freefree: both rotation blocks fully released → inactive (this is
        # the structural fix for the hinge/torsion singularity)
        frfr = N.FrameElement(n1, n2, sec; release=:freefree)
        sigf = N.dof_signature(frfr)
        @test all(sigf[1:3]) && all(sigf[7:9])
        @test !any(sigf[4:6]) && !any(sigf[10:12])

        # joist: torsion kept → whole rotation blocks stay active
        joist = N.FrameElement(n1, n2, sec; release=:joist)
        @test all(N.dof_signature(joist))

        # truss: translations only
        t = N.TrussElement(n1, n2, sec)
        sigt = N.dof_signature(t)
        @test sigt == (true, true, true, false, false, false,
            true, true, true, false, false, false)

        # signature is consistent with the stiffness matrix: unsigned slots
        # have identically zero rows/columns in GLOBAL coordinates
        for el in (ff, frfr, t)
            K = N.stiffness(el, n1.position, n2.position)
            s = N.dof_signature(el)
            for i in 1:12
                if !s[i]
                    @test all(iszero, K[i, :]) && all(iszero, K[:, i])
                end
            end
        end
    end

    @testset "stiffness delegates to kernels" begin
        el = N.FrameElement(n1, n2, sec)
        K = N.stiffness(el, n1.position, n2.position)
        @test Matrix(K) ≈ FIXTURES["fef_fixedfixed/global_K"] rtol = 1e-10

        t = N.TrussElement(n1, n2, sec)
        Kt = N.stiffness(t, n1.position, n2.position)
        B = N.truss_stiffness(sec, n1.position, n2.position)
        tidx = [1, 2, 3, 7, 8, 9]
        @test Matrix(Kt)[tidx, tidx] ≈ Matrix(B) rtol = 1e-12
    end

    @testset "element utilities" begin
        el = N.FrameElement(n1, n2, sec)
        p1, p2 = N.endpoints(el)
        @test p1 === n1.position && p2 === n2.position
        @test N.midpoint(el) ≈ (n1.position + n2.position) / 2
        Λ = N.local_frame(el)
        @test Matrix(Λ) ≈ FIXTURES["fef_fixedfixed/R"][1:3, 1:3] rtol = 1e-12
    end

    @testset "NodalSpring" begin
        base = N.Node([0.0, 0.0, 0.0], :free, :base)
        sp = N.NodalSpring(base, [0.0, 0.0, 5e4, 0.0, 0.0, 0.0], :soil)
        @test sp.node === base
        @test sp.stiffness[3] == 5e4

        uniform = N.NodalSpring(base, 1e5)
        @test collect(uniform.stiffness) == [1e5, 1e5, 1e5, 0.0, 0.0, 0.0]

        @test_throws AssertionError N.NodalSpring(base, [-1.0, 0, 0, 0, 0, 0])
        @test_throws AssertionError N.NodalSpring(base, [1.0, 0, 0])
    end

    @testset "show methods render" begin
        el = N.FrameElement(n1, n2, sec, :girder)
        plain = sprint(show, MIME"text/plain"(), el)
        @test occursin("girder", plain) && occursin("roll", plain) && occursin("EA", plain)

        t = N.TrussElement(n1, n2, sec)
        @test occursin("two-force", sprint(show, MIME"text/plain"(), t))

        sp = N.NodalSpring(N.Node([0.0, 0.0, 0.0], :free, :base), 1e5)
        plainsp = sprint(show, MIME"text/plain"(), sp)
        @test occursin("elastic support", plainsp) && occursin("k(Tx)", plainsp)
    end
end
