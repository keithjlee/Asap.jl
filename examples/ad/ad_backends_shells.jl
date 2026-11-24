#=
AD backends through shells — reproducible comparison.

Zygote, Enzyme (reverse) and Mooncake on the SAME shell walls: compliance
gradient w.r.t. all node positions, warm medians of 5. All three use Asap's
shell rules (local ForwardDiff Jacobians, `AsapShellADExt`); Mooncake and
Enzyme reach them through `AsapShellMooncakeExt` / `AsapShellEnzymeExt` —
both fail on the shell kernels natively. Enzyme needs runtime activity and
a `Const` function annotation (as everywhere in Asap).

Recorded 2026-11-24 (Julia 1.13.4, Enzyme 0.13.221, Mooncake 0.5.72, Apple M4):
    12 elements: primal 0.18 ms | Zygote 3.3 | Enzyme 3.1 | Mooncake 9.3 ms
    96 elements: primal 1.36 ms | Zygote 29.2 | Enzyme 25.6 | Mooncake 76.1 ms
(before the 1.3 linear-cost and split shell rules: Zygote 6.0 / 59.2, Enzyme
4.1 / 32.5, Mooncake 10.6 / 97.7 ms)
=#

begin
    using Asap, LinearAlgebra, Statistics
    using Zygote, ChainRulesCore, Enzyme, Mooncake, DifferentiationInterface
    import ForwardDiff
end

begin
    function wallmodel(size)
        conc = Material(3.0e4, 2.4e-3, 0.2)
        w = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 3.0, 2.0, ShellSection(conc, 0.2); size = size)
        top = [w.grid[i, end] for i in eachindex(w.s)]
        m = Model(w.nodes, AbstractElement{Float64}[w.elements...],
            AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81]) for el in w.elements]...,
                                  [NodeForce(n, [0.5, 0, 0]) for n in top]...])
        solve!(m)
        return m, extract_state(m), length(w.elements)
    end
end

begin
    rows = []                                    # collected for the report (last block)
    for sz in (0.75, 0.25)
        m, st, ne = wallmodel(sz)
        f(X) = compliance(m, ModelState{eltype(X)}(X, st.sections))
        tp = median([@elapsed f(st.X) for _ in 1:5])
        gz = Zygote.gradient(f, st.X)[1]
        tz = median([@elapsed Zygote.gradient(f, st.X) for _ in 1:5])
        be = AutoEnzyme(; mode = Enzyme.set_runtime_activity(Enzyme.Reverse), function_annotation = Enzyme.Const)
        ge = DifferentiationInterface.gradient(f, be, st.X)
        te = median([@elapsed DifferentiationInterface.gradient(f, be, st.X) for _ in 1:5])
        rule = Mooncake.build_rrule(f, st.X)
        gm = Mooncake.value_and_gradient!!(rule, f, st.X)[2][2]
        tm = median([@elapsed Mooncake.value_and_gradient!!(rule, f, st.X) for _ in 1:5])
        println("$ne elements: primal $(round(1e3tp, digits = 2)) ms | Zygote $(round(1e3tz, digits = 1)) ms | ",
            "Enzyme $(round(1e3te, digits = 1)) ms (rel. diff $(norm(ge - gz) / norm(gz))) | ",
            "Mooncake $(round(1e3tm, digits = 1)) ms (rel. diff $(norm(gm - gz) / norm(gz)))")
        push!(rows, (mesh_size = sz, elements = ne, nodes = length(m.nodes), n_vars = length(st.X),
            primal_ms = 1e3tp, zygote_ms = 1e3tz, enzyme_ms = 1e3te, mooncake_ms = 1e3tm,
            enzyme_rel_diff = norm(ge - gz) / norm(gz), mooncake_rel_diff = norm(gm - gz) / norm(gz)))
    end
end

#=
## Numbers for the report

`reports/D2_ad_backends_shells.typ` prints these. This script does not load
`figures.jl` (no plots), so it writes the JSON directly.
=#
begin
    using JSON
    _fin(x) = x isa AbstractFloat && !isfinite(x) ? nothing : x
    write(joinpath(@__DIR__, "..", "reports", "data", "D2.json"), JSON.json(Dict(
        "julia" => string(VERSION), "cpu" => Sys.CPU_NAME, "threads" => Threads.nthreads(),
        "versions" => Dict("enzyme" => string(pkgversion(Enzyme)), "zygote" => string(pkgversion(Zygote)),
            "mooncake" => string(pkgversion(Mooncake)), "differentiationinterface" => string(pkgversion(DifferentiationInterface))),
        "samples" => 5,
        "rows" => [Dict(string(k) => _fin(v) for (k, v) in pairs(r)) for r in rows]); pretty = true))
end
