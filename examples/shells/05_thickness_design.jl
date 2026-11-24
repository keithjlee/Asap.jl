#=
Walkthrough 5 — design sensitivities and thickness optimization

Shift+enter through the blocks. Units: kN, m.

Shells are differentiable end to end: wall thickness → stiffness AND self-
weight → displacements. With AsapOptim, `ThicknessVariable(elements, …)`
makes a whole wall (or a storey of it) one design variable.

The structure: a coupled wall line — a long pier (3 m) and a short pier
(1.5 m) joined by coupling beams — two storeys, lateral load at both
levels. Four variables: the thickness of each pier in each storey. We
compute the roof-drift gradient three ways (reverse AD, implicit
differentiation, finite differences), then redistribute concrete with a
classical optimality-criteria loop: minimum volume at the SAME drift.

Fast reverse AD through shells needs ForwardDiff loaded alongside Zygote
(Asap's shell rules); Mooncake works too (load it, same rules bridged).
=#

begin
    using Asap, AsapOptim, LinearAlgebra
    using Zygote, ForwardDiff
    include(joinpath(@__DIR__, "figures.jl"))
end

begin
    conc = Material(30.0e6, 2.4, 0.2)
    t0 = 0.25
    sec = ShellSection(conc, t0)
    H = 3.5
    levels = [H, 2H]
    long = mesh_wall([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 3.0, 2H, sec; z_lines = levels, size = 0.5, id = :long)
    short = mesh_wall([4.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.5, 2H, sec; z_lines = levels, size = 0.5, id = :short)
    nodeat(w, s, z) = only(filter(n -> n.position ≈ w.origin + s * w.along + z * w.up, w.nodes))
    cb = Section(conc, 0.12, 1.6e-3, 1.6e-3, 1.8e-3)
    cbeams = [FrameElement(nodeat(long, 3.0, z), nodeat(short, 0.0, z), cb, :coupling) for z in levels]
    roofnodes = vcat([long.grid[i, end] for i in eachindex(long.s)], [short.grid[i, end] for i in eachindex(short.s)])
    lvl1nodes = filter(n -> n.position[3] ≈ H, vcat(long.nodes, short.nodes))
    loads = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0.0, -9.81]) for el in vcat(long.elements, short.elements)]...,
        [NodeForce(n, [150.0 / length(lvl1nodes), 0, 0]) for n in lvl1nodes]...,
        [NodeForce(n, [250.0 / length(roofnodes), 0, 0]) for n in roofnodes]...]
    model = Model(vcat(long.nodes, short.nodes), AbstractElement{Float64}[long.elements..., short.elements..., cbeams...], loads)
    solve!(model)
    roof = long.grid[end, end]                       # roof node at the long pier's inner edge
    model
end

#=
## Design variables: one thickness per pier per storey
=#
begin
    storey(w, k) = [el for el in w.elements if levels[k] - H <= sum(n -> n.position[3], el.nodes) / 4 <= levels[k]]
    vars = AbstractVariable[ThicknessVariable(storey(long, 1), 0.12, 0.6), ThicknessVariable(storey(long, 2), 0.12, 0.6),
                            ThicknessVariable(storey(short, 1), 0.12, 0.6), ThicknessVariable(storey(short, 2), 0.12, 0.6)]
    names = ["long, storey 1", "long, storey 2", "short, storey 1", "short, storey 2"]
    p = OptParams(model, vars)
    idx = 6 * (roof.index - 1) + 1                   # roof ux in the full displacement vector
    drift(x) = solve_structure(x, p).U[idx]
    volume(x) = shell_volume(x, p)
    x0 = copy(p.values)
    (drift0_mm = 1000 * drift(x0), volume0 = volume(x0))
end

#=
## Gradients, three ways

- Zygote reverse mode through `solve_structure` (shell kernels use Asap's
  local-Jacobian rules),
- implicit differentiation: `solution_tangents` — one factorization and one
  multi-RHS back-substitution for ALL variables (the constraint-Jacobian
  workhorse),
- central finite differences, as the referee.
=#
begin
    gz = Zygote.gradient(drift, x0)[1]
    gi = solution_tangents(x0, p).dU[idx, :]
    h = 1e-6
    gfd = [(drift(x0 + h * (1:4 .== j)) - drift(x0 - h * (1:4 .== j))) / 2h for j in 1:4]
    grads = (zygote = gz, implicit = gi, finite_diff = gfd)
    @show grads
    f = fig_bars(names, -1000 .* gi; ylabel = "−∂drift/∂t [mm per m]",
        title = "where thickness buys stiffness")
    savefig("05_sensitivities.png", f)
end

#=
## Optimality criteria: minimum concrete at the starting drift

The classical resizing rule for one displacement constraint:
tᵢ ← tᵢ · (λ · (−∂δ/∂tᵢ) / (∂V/∂tᵢ))^½, with λ found by bisection so the
linearized drift equals the limit, and tᵢ clamped to its bounds. Each
iteration costs one implicit-tangent evaluation.
=#
begin
    δlim = drift(x0)
    lb, ub = p.lb, p.ub
    gV = Zygote.gradient(volume, x0)[1]              # ∂V/∂t = element areas (constant)
    function optimality_criteria(x0; iters = 25)
        x = copy(x0)
        hist = [(V = volume(x), δ = drift(x))]
        for it in 1:iters
            tg = solution_tangents(x, p)
            δ, g = tg.res.U[idx], tg.dU[idx, :]
            trial(λ) = clamp.(x .* (λ .* max.(-g, 1e-16) ./ gV) .^ 0.5, lb, ub)
            lin(λ) = δ + dot(g, trial(λ) - x)        # linearized drift
            lo, hi = 1e-12, 1e12
            for _ in 1:200
                mid = sqrt(lo * hi)
                lin(mid) > δlim ? (lo = mid) : (hi = mid)
            end
            x = trial(hi)
            push!(hist, (V = volume(x), δ = drift(x)))
        end
        return x, hist
    end
    x, hist = optimality_criteria(x0)
    result = (thickness = round.(x, digits = 3), volume_saving = 1 - hist[end].V / hist[1].V,
              drift_ratio_to_limit = hist[end].δ / δlim)
    @show result
end

begin
    f = Figure(size = (620, 380))
    ax1 = Axis(f[1, 1]; xlabel = "iteration", ylabel = "concrete volume [m³]", title = "optimality criteria history")
    ax2 = Axis(f[1, 1]; ylabel = "roof drift / limit", yaxisposition = :right)
    hidespines!(ax2); hidexdecorations!(ax2)
    l1 = scatterlines!(ax1, 0:length(hist)-1, [h.V for h in hist]; color = kjl_blue)
    l2 = scatterlines!(ax2, 0:length(hist)-1, [h.δ / δlim for h in hist]; color = kjl_pink)
    axislegend(ax1, [l1, l2], ["volume", "drift / limit"]; position = :rt)
    savefig("05_oc_history.png", f)
end

begin
    f = Figure(size = (620, 360))
    ax = Axis(f[1, 1]; ylabel = "thickness [m]", xticks = (1:4, names), xticklabelrotation = π / 8,
        title = "thickness: start vs optimized (same drift)")
    barplot!(ax, [1, 2, 3, 4] .- 0.18, x0; width = 0.34, color = kjl_gray, label = "start")
    barplot!(ax, [1, 2, 3, 4] .+ 0.18, x; width = 0.34, color = kjl_blue, label = "optimized")
    axislegend(ax; position = :rt)
    savefig("05_thickness.png", f)
end
