#=
Walkthrough O3 — buying exactly as much connection stiffness as you need

Shift+enter through the blocks. Units: kN, m.

A precast frame with PINNED column bases, so the beam-column connections
are the only lateral load path. Connection rotational stiffness costs money
(embedment, sleeves, welds). Design variables: the rotational stiffness of
every beam end (AsapOptim `JointVariable`). Problem: minimum total stiffness
(a proxy for connection cost) such that the roof drift ≤ H/400.

Each iteration needs the drift's gradient w.r.t. every joint: implicit
differentiation (`solution_tangents`) gives the whole solution tangent in
one factorization + one back-substitution.
=#

begin
    using Asap, AsapOptim, LinearAlgebra
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

#=
## A two-storey, three-bay precast frame
=#
begin
    conc = Material(30e6, 2.4, 0.2)
    col = Section(conc, 0.16, 2.13e-3, 2.13e-3, 3.6e-3)       # 400 × 400
    bm = Section(conc, 0.12, 1.6e-3, 0.9e-3, 2.0e-3)          # 300 × 400
    H, B, nb, ns = 3.5, 6.0, 3, 2
    grid = [Node([B * i, 0.0, H * j], j == 0 ? :pinned : :free) for i in 0:nb, j in 0:ns]
    cols = [FrameElement(grid[i+1, j+1], grid[i+1, j+2], col, :column) for i in 0:nb, j in 0:ns-1]
    beams = [FrameElement(grid[i+1, j+1], grid[i+2, j+1], bm, :beam) for i in 0:nb-1, j in 1:ns]
    loads = AbstractLoad{Float64}[[LineLoad(b, [0.0, 0.0, -25.0]) for b in beams]...,
        [NodeForce(grid[1, j+1], [20.0 * j / ns, 0.0, 0.0]) for j in 1:ns]...]
    frame = Model(vec(grid), AbstractElement{Float64}[vec(cols)..., vec(beams)...], loads)
    planarize!(frame, :ZX)
    solve!(frame)
    roof = grid[1, ns+1]
    limit = ns * H / 400
    (rigid_drift_mm = 1000 * displacement(frame.results, roof)[1], limit_mm = 1000 * limit)
end

#=
## Joint variables: both ends of every beam
=#
begin
    vars = AbstractVariable[JointVariable(b, pos, 2.0e5, 1.0e3, 1.0e8) for b in vec(beams) for pos in (:start, :end)]
    p = OptParams(frame, vars)
    idx = 6 * (roof.index - 1) + 1
    drift(x) = solve_structure(x, p).U[idx]
    (drift_at_start_mm = 1000 * drift(p.values), joints = length(vars))
end

#=
## Optimality criteria on log-stiffness

minimize Σ kᵢ  subject to  δ(k) ≤ δmax:
kᵢ ← kᵢ·(λ·(−∂δ/∂kᵢ))^½, λ by bisection so the linearized drift equals the
limit. −∂δ/∂kᵢ is where stiffness pays; the update moves stiffness there.
=#
begin
    function oc_joints(x0; iters = 40)
        x = copy(x0); hist = [(K = sum(x), δ = drift(x))]
        for _ in 1:iters
            t = solution_tangents(x, p)
            δ, g = t.res.U[idx], t.dU[idx, :]
            upd(λ) = clamp.(x .* (λ .* max.(-g, 1e-30)) .^ 0.5, p.lb, p.ub)
            lin(λ) = δ + dot(g, upd(λ) - x)
            lo, hi = 1e-30, 1e30
            for _ in 1:300
                mid = sqrt(lo * hi); lin(mid) > limit ? (lo = mid) : (hi = mid)
            end
            x = upd(hi); push!(hist, (K = sum(x), δ = drift(x)))
        end
        return x, hist
    end
    xk, hk = oc_joints(p.values)
    (total_stiffness_ratio = hk[end].K / hk[1].K, drift_over_limit = hk[end].δ / limit)
end

begin
    labels = ["beam $(i) $(pos == :start ? "start" : "end")" for i in 1:length(vec(beams)) for pos in (:start, :end)]
    f = fig_bars(labels, log10.(xk); ylabel = "log₁₀ k  [kN·m/rad]", rotation = π / 4,
        title = "optimized joint stiffness (drift = H/400): leeward beam ends carry it")
    savefig("O3_joints.png", f)
end

#=
## Numbers for the report

`reports/O3_frame_joints.pdf` prints these (see `figures.jl`, `report_data`).
=#
begin
    report_data("O3";
        E = conc.E, storey_height = H, bay = B, bays = nb, storeys = ns,
        col = (A = col.A, I = col.Ix), beam = (A = bm.A, I = bm.Ix),
        line_load = -loads[1].w[1] * loads[1].direction[3], lateral = [20.0 * j / ns for j in 1:ns],
        rigid_drift_mm = 1000 * displacement(frame.results, roof)[1], limit_mm = 1000 * limit,
        k0 = vars[1].value, k_lb = p.lb[1], k_ub = p.ub[1], joints = length(vars),
        drift_at_start_mm = 1000 * drift(p.values), iterations = length(hk) - 1,
        history_K = [h.K for h in hk], history_drift_mm = [1000 * h.δ for h in hk],
        total_stiffness_ratio = hk[end].K / hk[1].K, drift_over_limit = hk[end].δ / limit,
        joint = [(label = labels[k], bay = (k - 1) ÷ 2 % nb + 1, level = (k - 1) ÷ (2nb) + 1,
                  k = xk[k], at_lb = xk[k] <= 1.0001 * p.lb[k]) for k in eachindex(xk)])
end
