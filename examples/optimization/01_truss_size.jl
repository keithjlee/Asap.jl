#=
Walkthrough O1 — SIZE optimization of a truss (member areas)

Shift+enter through the blocks. Units: kN, m.

The Pratt truss of walkthrough A1, every member's area a design variable
(AsapOptim `AreaVariable`). Problem: minimum compliance (maximum stiffness)
for a fixed volume of steel. Gradients come from Zygote through the solve;
the optimizer is the classical optimality-criteria (OC) update, which is
exact-in-spirit for this problem: at the optimum every member has the same
strain-energy density — the truss becomes fully stressed.
=#

begin
    using Asap, AsapOptim, LinearAlgebra, Zygote
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
    include(joinpath(@__DIR__, "pratt.jl"))
end

begin
    truss, bottom, top, members = pratt_truss()
    vars = AbstractVariable[AreaVariable(el, 2.0e-3, 1.0e-6, 2.0e-2) for el in members]
    p = OptParams(truss, vars)
    x0 = copy(p.values)
    C(x) = compliance(solve_structure(x, p), p)
    L = [length(el) for el in members]
    V(x) = dot(x, L)
    V0 = V(x0)
    (compliance0 = C(x0), volume0 = V0)
end

#=
## Gradient check

Zygote's gradient against central finite differences (a few entries).
=#
begin
    g = Zygote.gradient(C, x0)[1]
    h = 1e-9
    fd = [(C(x0 + h * (eachindex(x0) .== j)) - C(x0 - h * (eachindex(x0) .== j))) / 2h for j in 1:4]
    (zygote = g[1:4], finite_diff = fd)
end

#=
## Optimality criteria

xᵢ ← clamp(xᵢ·(−∂C/∂xᵢ / (λ·Lᵢ))^½, lb, ub), with λ found by bisection so
the volume stays V0. (−∂C/∂Aᵢ = Nᵢ²Lᵢ/(E Aᵢ²) ≥ 0: compliance gradients of
truss areas are never positive.)
=#
begin
    function oc(x0; iters = 60)
        x = copy(x0); hist = Float64[C(x)]
        for _ in 1:iters
            g = Zygote.gradient(C, x)[1]
            upd(λ) = clamp.(x .* (max.(-g, 0.0) ./ (λ .* L)) .^ 0.5, p.lb, p.ub)
            lo, hi = 1e-20, 1e20
            for _ in 1:200
                mid = sqrt(lo * hi)
                V(upd(mid)) > V0 ? (lo = mid) : (hi = mid)
            end
            x = upd(hi); push!(hist, C(x))
        end
        return x, hist
    end
    xopt, hist = oc(x0)
    resopt = solve_structure(xopt, p)
    σ = abs.(axial_stress(resopt, p))
    active = xopt .> 10 * p.lb
    (compliance_reduction = 1 - hist[end] / hist[1], volume = V(xopt) / V0,
     stress_spread_active = extrema(σ[active]) ./ 1e3)            # [MPa]: ≈ uniform ⇒ fully stressed
end

begin
    f = fig_series(0:length(hist)-1, [("compliance", hist)]; xlabel = "OC iteration", ylabel = "compliance [kN·m]",
        title = "size optimization (volume fixed)", legendpos = :rt)
    savefig("O1_history.png", f)
end

begin
    opt = updatemodel(p, xopt)               # write the design into p.model (in place), re-solve, return it
    f = fig_lines(opt; plane = :XY, widths = 1 .+ 14 .* xopt ./ maximum(xopt),
        title = "size-optimized Pratt truss (width ∝ area; blue C / pink T)")
    savefig("O1_optimized.png", f)
end

#=
## Numbers for the report

`reports/O1_truss_size.pdf` prints these (see `figures.jl`, `report_data`).
=#
begin
    report_data("O1";
        span = bottom[end].position[1], depth = top[1].position[2], panels = length(bottom) - 1,
        n_members = length(members), panel_load = -truss.loads[1].value[2], n_loads = length(truss.loads),
        E = members[1].section.material.E, area0 = x0[1], lb = p.lb[1], ub = p.ub[1],
        compliance0 = hist[1], volume0 = V0, compliance_final = hist[end], volume_final = V(xopt),
        compliance_reduction = 1 - hist[end] / hist[1], volume_ratio = V(xopt) / V0,
        iterations = length(hist) - 1, history = hist,
        gradcheck = (h = h, zygote = g[1:4], fd = fd, rel_err = abs.(g[1:4] .- fd) ./ abs.(fd)),
        grad_max_positive = maximum(g),
        n_active = count(active), stress_active_MPa = collect(extrema(σ[active]) ./ 1e3),
        members = [(id = string(el.id), length = L[i], area_cm2 = 1e4 * xopt[i], stress_MPa = σ[i] / 1e3,
                    force = axial_force(resopt, p)[i], active = active[i]) for (i, el) in enumerate(members)])
end
