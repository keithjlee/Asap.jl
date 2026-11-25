#=
Walkthrough 4 — precast panels and their joints (LinkElement)

Shift+enter through the blocks. Units: kN, m.

A precast shear wall is not a monolith: storey-high panels sit on each
other at horizontal joints (grout bed + dowels or grouted sleeves), and
the joint's stiffness — not just the panels' — sets the wall's drift and
force flow. Here two 4 m × 3.5 m panels are stacked; the horizontal joint
is a line of zero-length `LinkElement`s between duplicate nodes (one on
each panel). We sweep the joint's shear stiffness and read the joint's
axial force distribution — where it needs tension ties.

This is LINEAR: a real joint that opens in tension is a contact/gap
problem (the nonlinear phase on Asap's roadmap). The linear result tells
you whether it would open, and how much tie force keeps it closed.
=#

begin
    using Asap, LinearAlgebra
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

begin
    conc = Material(30.0e6, 2.4, 0.2)
    panel = ShellSection(conc, 0.20)
    L, H = 4.0, 3.5
    Vtop = 100.0                                    # lateral force at the top [kN]
end

#=
## Two panels with duplicate nodes at the joint

The lower panel is fixed at its base; the upper panel is meshed with
`fixity_bottom = :free` so its bottom row is NEW nodes coincident with the
lower panel's top row. Links connect each coincident pair.
=#
begin
    function stacked(kjoint)
        lower = mesh_wall([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], L, H, panel; size = 0.25)
        upper = mesh_wall([0.0, 0.0, H], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], L, H, panel; size = 0.25,
            fixity_bottom = :free)
        joint = [LinkElement(lower.grid[i, end], upper.grid[i, 1], kjoint; id = :hjoint) for i in eachindex(lower.s)]
        top = [upper.grid[i, end] for i in eachindex(upper.s)]
        loads = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0.0, -9.81]) for el in vcat(lower.elements, upper.elements)]...,
            [NodeForce(n, [Vtop / length(top), 0.0, 0.0]) for n in top]...]
        m = Model(vcat(lower.nodes, upper.nodes), AbstractElement{Float64}[lower.elements..., upper.elements..., joint...], loads)
        solve!(m)
        drift = sum(n -> displacement(m.results, n)[1], top) / length(top)
        return (; m, lower, upper, joint, drift)
    end

    # link stiffness per joint node: (k along wall = shear, k through wall, k vertical = bearing,
    # kθ about the joint line, 0, 0). Every link sits ON the joint line, so
    # translational springs alone cannot stop the upper panel hinging
    # out-of-plane about it — a real joint resists that rotation through its
    # bearing width, modelled here by kθ.
    kn(kshear; kbear = 5e7, kθ = 1e5) = [kshear, kshear, kbear, kθ, 0.0, 0.0]
    stiff = stacked(kn(1e10; kbear = 1e10))         # effectively monolithic
    stiff.drift
end

#=
## Drift vs joint shear stiffness

Per-node shear stiffness from soft (a sliding joint) to rigid. The plateau
on the right is the monolithic wall; to the left, joint slip dominates.
=#
begin
    ks = 10.0 .^ (4:0.5:10)
    drifts = [stacked(kn(k)).drift for k in ks]
    f = fig_series(ks, [("two panels, linked", 1000 .* drifts)]; xscale = log10,
        xlabel = "joint shear stiffness per node [kN/m]", ylabel = "top drift [mm]",
        ref = 1000 * stiff.drift, reflabel = "monolithic", title = "drift vs joint stiffness",
        legendpos = :rt)
    savefig("04_joint_drift_sweep.pdf", f)
    @show extrema(drifts)
end

#=
## Force in the joint

With a realistic joint (shear 1e6 kN/m per node, bearing 5e7), the vertical
link forces are the joint's axial force distribution: compression (bearing)
at the toe, TENSION at the heel where the overturning exceeds the panel's
weight — that is the tie (dowel/sleeve) demand. The sum balances the
upper panel's weight exactly.
=#
begin
    r = stacked(kn(1e6))
    res = r.m.results
    # opening of each link (+ = tension): k·(u_upper − u_lower) along z
    s = [l.nodeStart.position[1] for l in r.joint]
    N = [l.k[3] * (displacement(res, l.nodeEnd)[3] - displacement(res, l.nodeStart)[3]) for l in r.joint]
    Wupper = areal_mass(panel) * 9.81 * L * H
    @show sum(N) (-Wupper)          # Σ joint forces = − weight of the upper panel
    f = Figure(size = (W_HALF, round(Int, 0.95W_HALF)))
    ax = Axis(f[1, 1]; xlabel = "position along joint s [m]", ylabel = "link axial force [kN]",
        title = "horizontal joint")
    barplot!(ax, s, N; color = force_color.(N), width = 0.2)
    hlines!(ax, [0.0]; color = kjl_black, linewidth = 0.6, linestyle = :solid)
    savefig("04_joint_forces.pdf", f)
end

begin
    f = fig_wall(r.m, r.upper; field = :Nyy, scale = 200.0, title = "upper panel: N_{yy}, ×200", width = W_HALF)
    savefig("04_upper_panel_Nyy.pdf", f)
end

#=
## Numbers for the report

`reports/S4_precast_joints.pdf` prints these (see `figures.jl`, `report_data`).
=#
begin
    report_data("S4";
        length = L, height = H, thickness = panel.t, E = conc.E, nu = conc.ν, rho = conc.ρ,
        Vtop = Vtop, mesh_size = 0.25, n_links = length(r.joint),
        n_elements_panel = length(r.lower.elements),
        k_bearing = 5e7, k_theta = 1e5, k_rigid = 1e10,
        monolithic_drift_mm = 1000 * stiff.drift,
        sweep = (k = collect(ks), drift_mm = 1000 .* drifts,
                 ratio = drifts ./ stiff.drift),
        realistic = (k_shear = 1e6, drift_mm = 1000 * r.drift, ratio = r.drift / stiff.drift),
        joint = (s = s, N = N, sum = sum(N), W_upper = Wupper,
                 max_tension = maximum(N), max_compression = minimum(N),
                 n_tension = count(>(0), N), tension_total = sum(filter(>(0), N)),
                 s_max_tension = s[argmax(N)], s_max_compression = s[argmin(N)]))
end
