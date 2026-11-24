#=
Walkthrough O2 — SHAPE vs SIZE optimization

Shift+enter through the blocks. Units: kN, m.

Same Pratt truss. SHAPE variables move the top-chord nodes vertically
(AsapOptim `SpatialVariable`, coupled symmetrically with `CoupledVariable`);
SIZE variables are the member areas (walkthrough O1). Which buys more?

A fair yardstick is C·V — compliance times volume. Scaling every area by s
scales V by s and C by 1/s, so C·V is independent of how much steel you
use: it measures the quality of the LAYOUT and the DISTRIBUTION of material.
=#

begin
    using Asap, AsapOptim, LinearAlgebra, Zygote
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
    include(joinpath(@__DIR__, "pratt.jl"))
end

#=
## Shape variables

Top nodes 1…5 (x = 2 … 10 m); symmetric pairs share a variable:
(1, 5), (2, 4) coupled, the midspan node 3 on its own. Each moves along Y
within ±1.5 m of its start.
=#
begin
    truss, bottom, top, members = pratt_truss()
    v1 = SpatialVariable(top[1], 0.0, -1.5, 1.5, :Y)
    v2 = SpatialVariable(top[2], 0.0, -1.5, 1.5, :Y)
    v3 = SpatialVariable(top[3], 0.0, -1.5, 1.5, :Y)
    shape_vars = AbstractVariable[v1, v2, v3, CoupledVariable(top[5], v1), CoupledVariable(top[4], v2)]
    ps = OptParams(truss, shape_vars)
    function CV(x)
        r = solve_structure(x, ps)
        return compliance(r, ps) * sum(r.A .* r.L)
    end
    CV0 = CV(ps.values)
end

#=
## Shape optimization: projected gradient descent with backtracking
=#
begin
    function descend(f, x0, lb, ub; iters = 80)
        x = copy(x0); hist = [f(x)]; step = 0.5
        for _ in 1:iters
            g = Zygote.gradient(f, x)[1]
            d = -g / norm(g)
            while step > 1e-6
                xn = clamp.(x .+ step .* d, lb, ub)
                if f(xn) < hist[end]
                    x = xn; push!(hist, f(xn)); step *= 1.5; break
                end
                step /= 2
            end
            step <= 1e-6 && break
        end
        return x, hist
    end
    xs, hs = descend(CV, ps.values, ps.lb, ps.ub)
    shaped = updatemodel(ps, xs)             # NB in place: `shaped === ps.model`, the truss itself
    (shape_only = hs[end] / CV0, top_chord_heights = [n.position[2] for n in shaped.nodes if n.id == :top])
end

#=
## Size optimization (from O1) on the ORIGINAL and the SHAPED truss

Min compliance at fixed volume = min C·V, so the same OC loop applies.
=#
begin
    function size_opt(model; iters = 60)
        els = model.elements
        p = OptParams(model, AbstractVariable[AreaVariable(el, el.section.A, 1e-6, 2e-2) for el in els])
        C(x) = compliance(solve_structure(x, p), p)
        L = [length(el) for el in els]; x = copy(p.values); V0 = dot(x, L)
        for _ in 1:iters
            g = Zygote.gradient(C, x)[1]
            upd(λ) = clamp.(x .* (max.(-g, 0.0) ./ (λ .* L)) .^ 0.5, p.lb, p.ub)
            lo, hi = 1e-20, 1e20
            for _ in 1:200
                mid = sqrt(lo * hi); dot(upd(mid), L) > V0 ? (lo = mid) : (hi = mid)
            end
            x = upd(hi)
        end
        return C(x) * dot(x, L), updatemodel(p, x), x
    end
    cv_size, sized, xsz = size_opt(pratt_truss()[1])
    cv_both, both, xb = size_opt(shaped)
    results = (initial = 1.0, size_only = cv_size / CV0, shape_only = hs[end] / CV0, shape_then_size = cv_both / CV0)
end

begin
    f = fig_bars(["initial", "size only", "shape only", "shape + size"], collect(values(results));
        ylabel = "C·V (normalized to initial)", title = "what shape and size each buy")
    savefig("O2_comparison.png", f)
end

begin
    f = fig_lines(both; plane = :XY, widths = 1 .+ 14 .* xb ./ maximum(xb),
        title = "shape + size optimized (width ∝ area)")
    savefig("O2_shape_size.png", f)
end

#=
## Numbers for the report

`reports/O2_truss_shape.typ` prints these (see `figures.jl`, `report_data`).
=#
begin
    report_data("O2";
        n_shape_vars = length(ps.values), shape_lb = ps.lb[1], shape_ub = ps.ub[1],
        CV0 = CV0, shape_x = xs, shape_history = hs ./ CV0, shape_iterations = length(hs) - 1,
        top_chord_heights = [n.position[2] for n in shaped.nodes if n.id == :top],
        top_chord_x = [n.position[1] for n in shaped.nodes if n.id == :top],
        top_chord_heights0 = [n.position[2] for n in pratt_truss()[3]],     # fresh: updatemodel moved `top`
        size_lb = 1e-6, size_ub = 2e-2, size_iterations = 60,
        cv_size = cv_size, cv_shape = hs[end], cv_both = cv_both,
        ratios = (initial = results.initial, size_only = results.size_only,
                  shape_only = results.shape_only, shape_then_size = results.shape_then_size),
        n_members = length(xb), n_at_lb_size = count(xsz .<= 10 * 1e-6), n_at_lb_both = count(xb .<= 10 * 1e-6))
end
