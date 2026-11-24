# Reverse rules for the frame kernels (ext/AsapFrameADExt.jl): the Zygote
# gradient through the plain-data kernels must equal ForwardDiff's, for every
# input class (rigidities, end springs incl. the exact 0/Inf limits,
# positions, roll angle), both section types and every element-load type —
# and the plain-data kernels must reproduce the struct kernels bitwise.

import ForwardDiff

@testset "frame kernel AD rules (AsapFrameADExt)" begin
    @test Base.get_extension(Asap, :AsapFrameADExt) !== nothing

    mat = Material(200e6, 77e6, 7.85, 0.3)
    secs = (Section(mat, 8e-3, 2e-4, 5e-5, 1e-6), RigiditySection(1.6e6, 4.0e4, 1.0e4, 77.0, 0.0628))
    endss = (EndConditions(:fixedfixed), EndConditions(:freefixed), EndConditions(:joist),
        EndConditions(EndSprings(Inf, Inf, 3e4, 2e4), EndSprings(5e6, Inf, Inf, 8e3)))
    x1, x2 = SVector(0.3, -0.2, 3.0), SVector(5.6, 0.4, 3.9)
    ψ = 0.7
    W = reshape(sin.(1:144), 12, 12)          # fixed cotangent
    # reference: ForwardDiff with the exact-limit springs (Inf rigid, 0
    # released) held constant, as the kernel's branches treat them — a dual
    # partial seeded on an Inf turns every direction into NaN, and a seeded 0
    # fails `iszero` and skips the release branch
    islimit(b) = !isfinite(b) || iszero(b)
    hold(kp, k) = map((a, b) -> islimit(b) ? oftype(a, b) : a, kp, k)

    @testset "stiffness: $(nameof(typeof(sec))), ends $i" for sec in secs, (i, ec) in enumerate(endss)
        r, k = Asap._stiffness_rigidities(sec), Asap._end_springs(ec)
        @test Asap._frame_K_plain(r, k, x1, x2, ψ) == Asap.frame_stiffness(sec, ec, x1, x2, ψ)
        f(r, k, a, b, θ) = sum(W .* Asap._frame_K_plain(r, k, a, b, θ))
        gz = Zygote.gradient(f, r, k, x1, x2, ψ)
        p0 = vcat(r, k, x1, x2, SVector(ψ))
        gf = ForwardDiff.gradient(p -> f(p[SVector(1, 2, 3, 4)], hold(p[SVector(5, 6, 7, 8, 9, 10, 11, 12)], k),
            p[SVector(13, 14, 15)], p[SVector(16, 17, 18)], p[19]), Vector(p0))
        @test vcat(gz[1], gz[2], gz[3], gz[4], SVector(gz[5])) ≈ gf rtol = 1e-12
        @test all(isfinite, vcat(gz[1], gz[2], gz[3], gz[4]))
    end

    @testset "element loads: $(nameof(typeof(sec)))" for sec in secs
        n1, n2 = Node(Vector(x1), :free), Node(Vector(x2), :free)
        el = FrameElement(n1, n2, sec)
        loads = (LineLoad(el, [0.0, 0.0, -10.0]), PointLoad(el, 0.3, [0.0, 4.0, -7.0]), SelfWeight(el))
        D = SVector{12}(cos.(1:12))
        for ld in loads, ec in endss, f in (Asap._fef_local_plain, Asap._fef_global_plain)
            r, k = Asap._fef_rigidities(sec), Asap._end_springs(ec)
            h(r, k, a, b, θ) = sum(D .* f(ld, r, k, a, b, θ))
            gz = Zygote.gradient(h, r, k, x1, x2, ψ)
            p0 = vcat(r, k, x1, x2, SVector(ψ))
            gf = ForwardDiff.gradient(p -> h(p[SVector(1, 2, 3)], hold(p[SVector(4, 5, 6, 7, 8, 9, 10, 11)], k),
                p[SVector(12, 13, 14)], p[SVector(15, 16, 17)], p[18]), Vector(p0))
            @test vcat(gz[1], gz[2], gz[3], gz[4], SVector(gz[5])) ≈ gf rtol = 1e-12 atol = 1e-12
        end
        ld = loads[1]
        @test Asap.element_fef_global(ld, sec, endss[4], x1, x2, ψ) ==
              Asap._fef_global_plain(ld, Asap._fef_rigidities(sec), Asap._end_springs(endss[4]), x1, x2, ψ)
    end

    @testset "model gradient through the rules matches ForwardDiff" begin
        # portal frame, line load + self-weight on the beam, semi-rigid beam ends
        sec = secs[1]
        n = [Node([0.0, 0, 0], :fixed), Node([0.0, 0, 3], :free), Node([5.0, 0, 3.2], :free), Node([5.0, 0, 0], :fixed)]
        beam = FrameElement(n[2], n[3], sec, EndConditions(EndSprings(Inf, Inf, Inf, 4e4), EndSprings(Inf, Inf, Inf, 2e4)), :beam)
        m = Model(n, AbstractElement{Float64}[FrameElement(n[1], n[2], sec), beam, FrameElement(n[4], n[3], sec)],
            AbstractLoad{Float64}[LineLoad(beam, [0.0, 0.0, -10.0]), SelfWeight(beam), NodeForce(n[2], [5.0, 0, 0])])
        solve!(m)
        st = extract_state(m)
        f(X) = compliance(m, ModelState{eltype(X)}(X, st.sections))
        @test Zygote.gradient(f, st.X)[1] ≈ ForwardDiff.gradient(f, st.X) rtol = 1e-10
    end
end
