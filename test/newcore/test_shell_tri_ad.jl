# S7g: triangles on the differentiable path — the ForwardDiff kernel rules
# (AsapShellADExt) against direct ForwardDiff, and model gradients (positions,
# thickness, design-dependent triangle loads) through the linear-cost group rules
# (AsapChainRulesExt) against finite differences.

using LinearAlgebra, StaticArrays, Zygote, FiniteDifferences, ChainRulesCore
import ForwardDiff

@testset "Triangles on the differentiable path (S7g)" begin
    sec = ShellSection(Material(3.0e4, 2.4e-3, 0.2), 0.2)
    A, D, S, γ = Asap._shell_rigidities(sec)
    X = (SVector(0.1, -0.2, 0.3), SVector(1.9, 0.1, 0.5), SVector(0.7, 1.4, -0.2))

    @testset "kernel rules ≡ ForwardDiff" begin
        W = reshape(sin.(1:324), 18, 18)
        pk(v) = (SMatrix{3,3}(v[1:9]), SMatrix{3,3}(v[10:18]), SMatrix{2,2}(v[19:22]), v[23],
                 SVector{3}(v[24:26]), SVector{3}(v[27:29]), SVector{3}(v[30:32]))
        h(args...) = sum(W .* Asap._shell_K_plain(args...))
        gz = Zygote.gradient(h, A, D, S, γ, X...)
        gf = ForwardDiff.gradient(v -> h(pk(v)...), Vector(vcat(vec(A), vec(D), vec(S), γ, X...)))
        @test vcat(vec(gz[1]), vec(gz[2]), vec(gz[3]), gz[4], gz[5], gz[6], gz[7]) ≈ gf rtol = 1e-12
        @test gz[4] == 0                                           # ANDES: no drilling penalty
        d = cos.(1:18)
        qg, ql = SVector(0.0, 0.0, -2.0), SVector(0.1, 0.0, 0.3)
        hf(a, b, c, X1, X2, X3) = dot(d, Asap._shell_fef_plain(a, b, c, X1, X2, X3))
        gz = Zygote.gradient(hf, A, qg, ql, X...)
        pf(v) = (A, SVector{3}(v[1:3]), SVector{3}(v[4:6]), SVector{3}(v[7:9]), SVector{3}(v[10:12]), SVector{3}(v[13:15]))
        gf = ForwardDiff.gradient(v -> hf(pf(v)...), Vector(vcat(qg, ql, X...)))
        @test vcat(gz[2], gz[3], gz[4], gz[5], gz[6]) ≈ gf rtol = 1e-12
    end

    @testset "model gradients: positions, thickness, triangle surface loads" begin
        w = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 3.0, 2.0, sec; size = 0.5)
        tris = AbstractElement{Float64}[t for el in w.elements for t in
            (ShellElement((el.nodes[1], el.nodes[2], el.nodes[3]), sec),
             ShellElement((el.nodes[1], el.nodes[3], el.nodes[4]), sec))]
        top = [w.grid[i, end] for i in eachindex(w.s)]
        loads = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81]) for el in tris]...,
            [AreaLoad(el, 0.05) for el in tris[1:3:end]]...,
            [NodeForce(n, [1.0, 0.2, 0]) for n in top]...]
        m = Model(w.nodes, tris, loads)
        solve!(m)
        st = extract_state(m)
        cfg = Zygote.ZygoteRuleConfig()
        g = only(m.cache.groups)
        @test rrule(cfg, Asap._group_entries_nnode_all, g, st, Val(3)) !== nothing
        @test rrule(cfg, Asap._shell3_fef, m.cache.fef, st) !== nothing
        nX = length(st.X)
        idx = eachindex(tris)
        function f(x)
            X = reshape(x[1:nX], 3, :)
            secs = map(i -> ShellSection(sec.material, x[nX+i]), idx)
            return compliance(m, ModelState{Float64}(X, secs))
        end
        x0 = vcat(vec(st.X), 0.2 .+ 0.01 .* sin.(idx))
        gz = Zygote.gradient(f, x0)[1]
        gfd = FiniteDifferences.grad(central_fdm(5, 1), f, x0)[1]
        @test gz ≈ gfd rtol = 1e-6
        @test norm(gz[nX+1:end]) > 0
    end
end
