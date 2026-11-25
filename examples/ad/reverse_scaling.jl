#=
Reverse-mode cost against model size — frames and shells.

Zygote compliance gradients with respect to every node coordinate, on grid
frames with line loads on every beam and on shell walls with self-weight,
next to the primal solve and a ForwardDiff check.

Recorded 2026-11-24 (Julia 1.13.4, Apple Silicon, one thread). "v1.2" is the
code path before any 1.3 rule; "1.3 kernels" adds the local-ForwardDiff kernel
rules only (AsapFrameADExt, the first AsapShellADExt); "1.3" is the release:
linear-cost group rules (AsapChainRulesExt) and the split shell rules.

    model                       primal   v1.2       1.3 kernels   1.3
    frame 120 el / 64 nodes     0.6 ms    62.8 ms     16.7 ms      7.1 ms
    frame 665 el / 294 nodes    2.9 ms   435.6 ms    256.2 ms     37.1 ms
    shell 96 el                 1.3 ms      —         63.9 ms     28.9 ms
    shell 384 el                5.5 ms      —        294.7 ms    104.4 ms
    shell 1536 el              24.0 ms      —       1976.8 ms    486.1 ms

With the group rules the cost is a FIXED multiple of the primal (frames 13×,
shells about 20×): the per-element closures index arrays they capture, and under
plain Zygote each index pulled back into a full-size cotangent (quadratic
cost). The rules leave the primal untouched (Enzyme and Mooncake never see
them) and accumulate one cotangent per state field.
=#

begin
    using Asap, LinearAlgebra, Statistics
    using Zygote, ChainRulesCore
    import ForwardDiff
end

begin
    function shellwall(sz)
        conc = Material(3.0e4, 2.4e-3, 0.2)
        w = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 6.0, 4.0, ShellSection(conc, 0.2); size = sz)
        top = [w.grid[i, end] for i in eachindex(w.s)]
        m = Model(w.nodes, AbstractElement{Float64}[w.elements...],
            AbstractLoad{Float64}[[NodeForce(n, [0.5, 0, 0]) for n in top]...,
                                  [SelfWeight(el; g = [0.0, 0, -9.81]) for el in w.elements]...])
        solve!(m)
        return m
    end

    # nx × ny bays of 6 m, nz storeys of 4 m, 10 kN/m on every beam, wind at the roof
    function gridframe(nx, ny, nz)
        mat = Material(200e6, 77e6, 7.85, 0.3)
        col = Section(mat, 1e-2, 1e-4, 1e-4, 2e-6)
        bm = Section(mat, 8e-3, 2e-4, 5e-5, 1e-6)
        node = Dict{NTuple{3,Int},Node{Float64}}()
        getn(i, j, k) = get!(() -> Node([6.0i, 6.0j, 4.0k], k == 0 ? :fixed : :free), node, (i, j, k))
        els = AbstractElement{Float64}[]
        loads = AbstractLoad{Float64}[]
        for i in 0:nx, j in 0:ny, k in 1:nz
            push!(els, FrameElement(getn(i, j, k - 1), getn(i, j, k), col))
            for (di, dj) in ((1, 0), (0, 1))
                (i + di <= nx && j + dj <= ny) || continue
                b = FrameElement(getn(i, j, k), getn(i + di, j + dj, k), bm)
                push!(els, b)
                push!(loads, LineLoad(b, [0.0, 0.0, -10.0]))
            end
        end
        append!(loads, [NodeForce(getn(i, 0, nz), [5.0, 0, 0]) for i in 0:nx])
        m = Model(collect(values(node)), els, loads)
        solve!(m)
        return m
    end
end

begin
    rows = []
    for (name, m) in (("frame 3×3×3", gridframe(3, 3, 3)), ("frame 6×6×5", gridframe(6, 6, 5)),
                      ("shell 0.5 m", shellwall(0.5)), ("shell 0.25 m", shellwall(0.25)), ("shell 0.125 m", shellwall(0.125)))
        st = extract_state(m)
        f(X) = compliance(m, ModelState{eltype(X)}(X, st.sections))
        g = Zygote.gradient(f, st.X)[1]
        tz = median([@elapsed Zygote.gradient(f, st.X) for _ in 1:3])
        tp = median([@elapsed f(st.X) for _ in 1:5])
        check = length(st.X) < 1000 ? norm(g - ForwardDiff.gradient(f, st.X)) / norm(g) : NaN
        push!(rows, (; name, elements = length(m.elements), nodes = length(m.nodes), primal_ms = 1e3tp, zygote_ms = 1e3tz, check))
        println(rpad(name, 14), " $(length(m.elements)) el: primal $(round(1e3tp, digits = 2)) ms, ",
            "Zygote $(round(1e3tz, digits = 1)) ms ($(round(tz / tp, digits = 0))× primal)",
            isnan(check) ? "" : ", vs ForwardDiff $(round(check, sigdigits = 2))")
    end
end
