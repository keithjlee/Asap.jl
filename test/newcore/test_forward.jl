# Forward-mode differentiation: the AsapForwardDiffExt Dual solve path, the
# solve_free frule, and the constructor promotion that lets Dual numbers
# traverse sections. Reverse-mode (Zygote) is the reference throughout —
# it is itself FD-verified in test_ad.jl.

using ForwardDiff
using ForwardDiff: Dual, value, partials
import ChainRulesCore

@testset "forward mode" begin

    @testset "Section/Material promotion" begin
        mat = Material(200e6, 1.0, 80.0, 0.3)
        d = Dual{:tag}(1e-2, 1.0)
        s = Section(mat, d, 1e-5, 2e-5, 3e-6)
        @test s isa Section{typeof(d)}
        @test s.material isa Material{typeof(d)}           # material lifted too
        @test value(EA(s)) ≈ 200e6 * 1e-2
        @test partials(EA(s), 1) ≈ 200e6                   # ∂(E·A)/∂A = E
        # Float64 stays Float64 (no promotion tax on the normal path)
        @test Section(mat, 1e-2) isa Section{Float64}
        @test Section(mat, 1e-2).material === mat
    end

    @testset "Dual solve_free ≡ implicit differentiation" begin
        K0 = sparse([1, 2, 2, 1, 3, 3], [1, 1, 2, 2, 3, 2],
            [4.0, 1.0, 3.0, 1.0, 2.0, 0.5], 3, 3)
        K0 = sparse(K0 + K0') / 2 + 3I                      # SPD
        F = [1.0, -2.0, 0.5]

        # seed ∂/∂nz₁: K(θ) = K0 + θ·e₁e₁ᵀ-ish through the nzval
        seed = zeros(nnz(K0)); seed[1] = 1.0
        nz_d = Dual{:t}.(K0.nzval, seed)
        Kd = SparseMatrixCSC(3, 3, K0.colptr, K0.rowval, nz_d)

        u = K0 \ collect(F)
        du_analytic = -(K0 \ (sparse(K0.rowval[1:1], [1], [1.0], 3, 3) * u))

        ud = Asap.solve_free(Kd, F)                        # Dual K, plain F
        @test value.(ud) ≈ u rtol = 1e-12
        @test [partials(x, 1) for x in ud] ≈ du_analytic rtol = 1e-12

        # multi-RHS (matrix) path
        Fm = [F 2 .* F]
        um = Asap.solve_free(Kd, Fm)
        @test value.(um) ≈ K0 \ Fm rtol = 1e-12
        @test [partials(x, 1) for x in um[:, 2]] ≈ 2 .* du_analytic rtol = 1e-12

        # plain K, Dual F: u̇ = K⁻¹Ḟ
        Fd = Dual{:t}.(F, ones(3))
        uf = Asap.solve_free(K0, Fd)
        @test [partials(x, 1) for x in uf] ≈ K0 \ ones(3) rtol = 1e-12

        # Dual K AND Dual F
        ub = Asap.solve_free(Kd, Fd)
        @test [partials(x, 1) for x in ub] ≈ du_analytic .+ K0 \ ones(3) rtol = 1e-12
    end

    @testset "solve_free frule ≡ rrule (transpose check)" begin
        K0 = sparse([1, 2, 2, 1], [1, 1, 2, 2], [4.0, 1.0, 1.0, 3.0])
        F = [1.0, 2.0]
        ΔK = SparseMatrixCSC(2, 2, K0.colptr, K0.rowval, [0.3, -0.1, -0.1, 0.2])
        ΔF = [0.5, -1.0]
        ū = [2.0, -0.7]

        u1, u̇ = ChainRulesCore.frule((ChainRulesCore.NoTangent(), ΔK, ΔF),
            Asap.solve_free, K0, F)
        u2, pb = ChainRulesCore.rrule(Asap.solve_free, K0, F)
        _, K̄, F̄ = pb(ū)
        @test u1 ≈ u2
        # ⟨ū, u̇⟩ = ⟨K̄, ΔK⟩ + ⟨F̄, ΔF⟩ — the defining adjoint identity
        @test dot(ū, u̇) ≈ dot(ChainRulesCore.unthunk(K̄), ΔK) + dot(ChainRulesCore.unthunk(F̄), ΔF) rtol = 1e-12
    end

    @testset "ForwardDiff end-to-end: pure solve vs Zygote" begin
        mat = Material(200e6, 1.0, 80.0, 0.3)
        sec = Section(mat, 1e-2)
        rot = [true, true, true]
        n1 = Node([0.0, 0.0, 0.0], vcat([false, false, false], rot))
        n2 = Node([4.0, 0.0, 0.0], vcat([true, true, false], rot))
        n3 = Node([8.0, 3.0, 0.0], vcat([false, false, false], rot))
        els = AbstractElement{Float64}[
            TrussElement(n1, n2, sec), TrussElement(n2, n3, sec)]
        model = Model([n1, n2, n3], els,
            AbstractLoad{Float64}[NodeForce(n2, [0.0, -50.0, 0.0])])
        solve!(model)
        st = Asap.extract_state(model)

        # θ: [Δy of node 2, area scale] — geometry AND section derivatives
        δ = zeros(size(st.X)); δ[2, 2] = 1.0               # constant seed matrix
        function obj(θ)
            X = st.X + θ[1] * δ                            # pure (Zygote-safe)
            sections = [Section(mat, 1e-2 * θ[2]), Section(mat, 1e-2 * θ[2])]
            ea = [EA(s) for s in sections]
            state = Asap.ModelState{eltype(θ)}(X, sections, ea, nothing)
            u = Asap.solve(model, state)
            return dot(u, u)
        end
        θ0 = [0.2, 1.1]
        g_fwd = ForwardDiff.gradient(obj, θ0)
        g_rev = Zygote.gradient(obj, θ0)[1]
        @test g_fwd ≈ g_rev rtol = 1e-10
    end
end
