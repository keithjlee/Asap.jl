# Linear-cost reverse rules on the per-element maps (ext/AsapChainRulesExt.jl):
# frame and shell stiffness groups, design-dependent frame and shell loads,
# `fef_local`. Under Zygote the pure path's per-element closures would pull
# every captured-array index back into a full-size cotangent (O(n²)); the rules
# differentiate each element on its own arguments and accumulate one cotangent
# per state field. Checked here against ForwardDiff (frames) and finite
# differences (shells) for every state field the rules touch: positions,
# sections (both shell section types), connection springs; several loads on
# one member; inactive slots; every shell surface-load type.

using Zygote, FiniteDifferences, StaticArrays, LinearAlgebra, ChainRulesCore
import ForwardDiff

@testset "group AD rules (linear-cost reverse)" begin
    cfg = Zygote.ZygoteRuleConfig()

    @testset "frames: positions, sections, end springs, several loads, fef_local" begin
        mat = Material(200e6, 77e6, 7.85, 0.3)
        sec = Section(mat, 8e-3, 2e-4, 5e-5, 1e-6)
        n = [Node([0.0, 0, 0], :fixed), Node([0.0, 0, 3], :free), Node([5.0, 0, 3.2], :free),
             Node([5.0, 0, 0], :fixed), Node([5.0, 4.0, 3.1], :free), Node([5.0, 4.0, 0], :fixed),
             Node([8.0, 4.0, 3.0], :free)]
        semi = EndConditions(EndSprings(Inf, Inf, Inf, 4e4), EndSprings(Inf, Inf, 3e4, 2e4))
        b1 = FrameElement(n[2], n[3], sec, semi, :beam)
        b2 = FrameElement(n[3], n[5], sec, EndConditions(:freefixed), :beam)
        b3 = FrameElement(n[5], n[7], sec, EndConditions(:fixedfree), :beam)   # n7 rotations: inactive
        els = AbstractElement{Float64}[FrameElement(n[1], n[2], sec), b1, FrameElement(n[4], n[3], sec),
            b2, FrameElement(n[6], n[5], sec), b3]
        loads = AbstractLoad{Float64}[LineLoad(b1, [0.0, 0.0, -10.0]), PointLoad(b1, 0.3, [0.0, 4.0, -7.0]),
            SelfWeight(b1), LineLoad(b2, [1.0, 0.0, -6.0]), PointLoad(b3, 0.6, [0.0, 0.0, -5.0]),
            NodeForce(n[2], [5.0, 0, 0])]
        m = Model(n, els, loads)
        solve!(m)
        st = extract_state(m)

        g = only(filter(g -> eltype(g.elements) <: FrameElement, m.cache.groups))
        @test any(s -> length(s) < 12, g.slots)                       # inactive slots exercised
        @test rrule(cfg, Asap._group_entries_all, g, st) !== nothing   # the rule is what Zygote uses
        @test rrule(cfg, Asap._frame_fef, m.cache.fef, st) !== nothing

        W = reshape(sin.(1:12length(m.cache.fef.loads)), 12, :)
        # x = [positions; areas; the four finite connection springs of b1]
        nX = length(st.X)
        x0 = vcat(vec(st.X), fill(8e-3, length(els)), [4e4, 3e4, 2e4])
        # the ends of every other element are the model's: rebuild them exactly
        ends0 = [g.endss[findfirst(==(i), g.model_indices)] for i in eachindex(els)]
        function obj(x)
            T = eltype(x)
            X = reshape(x[1:nX], 3, :)
            secs = map(i -> Section(mat, x[nX+i], 2e-4, 5e-5, 1e-6), eachindex(els))
            k = x[nX+length(els)+1:end]
            ends = map(i -> i == 2 ? EndConditions(EndSprings(Inf, Inf, Inf, k[1]), EndSprings(Inf, Inf, k[2], k[3])) :
                       ends0[i], eachindex(els))
            s = ModelState{T}(X, secs, nothing, ends)
            return compliance(m, s) + 1e-3 * sum(W .* fef_local(m.cache, s))
        end
        @test obj(x0) ≈ compliance(m, st) + 1e-3 * sum(W .* fef_local(m.cache, st)) rtol = 1e-12
        gz = Zygote.gradient(obj, x0)[1]
        gf = ForwardDiff.gradient(obj, x0)
        @test gz ≈ gf rtol = 1e-10
        @test all(isfinite, gz)
    end

    @testset "shells: positions, ShellSection + ShellRigiditySection, all surface loads" begin
        mat = Material(3.0e4, 2.4e-3, 0.2)
        sec = ShellSection(mat, 0.2)
        w = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 2.0, 1.5, sec; size = 0.5)
        nel = length(w.elements)
        loads = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81]) for el in w.elements]...,
            [AreaLoad(el, 0.05) for el in w.elements[1:2:end]]...,                         # local pressure
            [AreaLoad(el, [0.02, 0.0, -0.03]) for el in w.elements[2:2:end]]...]           # global traction
        top = [w.grid[i, end] for i in eachindex(w.s)]
        append!(loads, [NodeForce(nd, [1.0, 0.3, 0]) for nd in top])
        m = Model(w.nodes, AbstractElement{Float64}[w.elements...], loads)
        solve!(m)
        st = extract_state(m)
        g = only(m.cache.groups)
        @test rrule(cfg, Asap._group_entries_nnode_all, g, st, Val(4)) !== nothing
        @test rrule(cfg, Asap._shell_fef, m.cache.fef, st) !== nothing

        nX = length(st.X)
        # odd elements: ShellSection(t); even: cracked ShellRigiditySection of ShellSection(t)
        mksec(i, t) = isodd(i) ? ShellSection(mat, t) : ShellRigiditySection(ShellSection(mat, t); m11 = 0.35, m22 = 0.5, m12 = 0.35)
        function obj(x)
            X = reshape(x[1:nX], 3, :)
            secs = map(i -> mksec(i, x[nX+i]), 1:nel)
            return compliance(m, ModelState{Float64}(X, secs))
        end
        x0 = vcat(vec(st.X), 0.2 .+ 0.01 .* sin.(1:nel))
        gz = Zygote.gradient(obj, x0)[1]
        gfd = FiniteDifferences.grad(central_fdm(5, 1), obj, x0)[1]
        @test gz ≈ gfd rtol = 1e-6
        @test norm(gz[nX+1:end]) > 0                       # thickness reaches both section types

        # ShellSection → rigidities rule against ForwardDiff through the accessors
        Δ = (SMatrix{3,3}(sin.(1:9)), SMatrix{3,3}(cos.(1:9)), SMatrix{2,2}(1.0, -2.0, 0.5, 3.0), 0.7)
        h(p) = let s = ShellSection(Material(p[1], p[2], mat.ρ, p[3]), p[4]), r = Asap._shell_rigidities(s)
            sum(Δ[1] .* r[1]) + sum(Δ[2] .* r[2]) + sum(Δ[3] .* r[3]) + Δ[4] * r[4]
        end
        p0 = [mat.E, mat.G, mat.ν, 0.2]
        _, pb = rrule(Asap._shell_rigidities, ShellSection(mat, 0.2))
        ts = pb(Δ)[2]
        @test [ts.material.E, ts.material.G, ts.material.ν, ts.t] ≈ ForwardDiff.gradient(h, p0) rtol = 1e-12
    end
end
