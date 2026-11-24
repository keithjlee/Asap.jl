#=
Walkthrough A5 — semi-rigid joints, spring supports, variable sections

Shift+enter through the blocks. Units: kN, m.

Three ways real structures differ from "rigid or pinned, fixed or free,
prismatic":
1. semi-rigid beam-column connections (`EndSprings` — finite rotational
   stiffness between the classical limits),
2. elastic supports (`NodalSpring` — a footing on soil),
3. members whose section varies along their length (`VariableElement` — a
   haunched girder).
=#

begin
    using Asap, LinearAlgebra
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

begin
    steel = Material(200e6, 77e6, 7.85, 0.3)
    column = Section(steel, 1.3e-2, 3.0e-4, 1.0e-4, 1.0e-6)
    girder = Section(steel, 1.1e-2, 4.5e-4, 1.2e-4, 1.0e-6)
end

#=
## 1. Semi-rigid connections

`EndSprings(kx, kt, ky, kz)` per end (axial, torsion, rotations about the
local y and z axes; `Inf` = rigid, `0` = released). We sweep the girder's
rotational connection stiffness k from pinned to rigid and watch the
portal's drift and the girder's end moment. The classical releases are the
exact limits.
=#
begin
    function portal(ends)
        a = Node([0.0, 0, 0], :pinned); b = Node([0.0, 0, 4.0], :free)
        c = Node([8.0, 0, 4.0], :free); d = Node([8.0, 0, 0], :pinned)
        beam = FrameElement(b, c, girder, ends, :girder)
        m = Model([a, b, c, d], AbstractElement{Float64}[FrameElement(a, b, column), beam,
                FrameElement(d, c, column)],
            AbstractLoad{Float64}[LineLoad(beam, [0.0, 0.0, -20.0]), NodeForce(b, [15.0, 0, 0])])
        planarize!(m, :ZX)
        solve!(m)
        return 1000 * displacement(m.results, b)[1]
    end
    spring(k) = EndConditions(EndSprings(Inf, Inf, k, k), EndSprings(Inf, Inf, k, k))
    ks = 10.0 .^ (3:0.5:8)
    drifts = [portal(spring(k)) for k in ks]
    rigid = portal(EndConditions(:fixedfixed))
    f = fig_series(ks, [("semi-rigid girder joints", drifts)]; xscale = log10, legendpos = :rt,
        xlabel = "connection stiffness k [kN·m/rad]", ylabel = "drift [mm]", ref = rigid, reflabel = "rigid limit",
        title = "pinned-base portal: drift vs joint stiffness")
    savefig("A5_joint_sweep.png", f)
    (k_1e3 = drifts[1], k_1e8 = drifts[end], rigid = rigid)
end

#=
## 2. Spring supports

A `NodalSpring` is an elastic support attached to a node (it lives on the
model, like a load). A column on a footing: translational springs from the
soil, and a rotational spring kθ from the footing's rocking stiffness.
=#
begin
    function cantilever(kθ)
        base = Node([0.0, 0, 0], :free); top = Node([0.0, 0, 4.0], :free)
        el = FrameElement(base, top, column)
        m = Model([base, top], AbstractElement{Float64}[el], AbstractLoad{Float64}[NodeForce(top, [10.0, 0, -50.0])];
            springs = [NodalSpring(base, [1e6, 1e6, 1e6, kθ, kθ, 1e6], :footing)])
        planarize!(m, :ZX)
        solve!(m)
        return 1000 * displacement(m.results, top)[1]
    end
    kθs = 10.0 .^ (3:0.5:8)
    EI = 200e6 * 3.0e-4
    # cantilever bending + rigid rotation at the base + slide on the soil spring
    exact(kθ) = 1000 * (10 * 4^3 / (3EI) + 10 * 4 * 4 / kθ + 10 / 1e6)
    tops = cantilever.(kθs)
    f = fig_series(kθs, [("Asap", tops), ("closed form", exact.(kθs))]; xscale = log10, legendpos = :rt,
        xlabel = "footing rotational stiffness kθ [kN·m/rad]", ylabel = "top displacement [mm]",
        title = "column on a rotational spring")
    savefig("A5_footing_spring.png", f)
    maximum(abs.(tops .- exact.(kθs)) ./ exact.(kθs))     # ≈ 0: Asap matches the closed form
end

#=
## 3. A haunched girder (VariableElement)

One member, three segments: deep haunches over the supports (0–20% and
80–100% of the span) and a shallower midspan. Interior joints become
internal DOFs — no extra nodes in your model — and you query the member as
one piece.
=#
begin
    deep = Section(steel, 1.6e-2, 9.0e-4, 2.0e-4, 2.0e-6)
    a = Node([0.0, 0, 0], :fixed); b = Node([12.0, 0, 0], :fixed)
    haunched = VariableElement(a, b, AbstractSection{Float64}[deep, girder, deep], [0.2, 0.8], :haunched)
    hm = Model([a, b], AbstractElement{Float64}[haunched], AbstractLoad{Float64}[LineLoad(haunched, [0.0, 0.0, -25.0])])
    planarize!(hm, :ZX)
    solve!(hm)
    prismatic = FrameElement(Node([0.0, 0, 0], :fixed), Node([12.0, 0, 0], :fixed), girder)
    pm = Model([prismatic.nodeStart, prismatic.nodeEnd], AbstractElement{Float64}[prismatic],
        AbstractLoad{Float64}[LineLoad(prismatic, [0.0, 0.0, -25.0])])
    planarize!(pm, :ZX); solve!(pm)
    pick(f) = maximum(abs, f.Mz) > maximum(abs, f.My) ? f.Mz : f.My
    Mh = pick(InternalForces(hm, haunched; resolution = 60)); Mp = pick(InternalForces(pm, prismatic; resolution = 60))
    # stiffer ends attract moment: larger support moments, smaller midspan moment than wL²/24
    (haunched_support = maximum(abs, Mh[[1, end]]), prismatic_support = maximum(abs, Mp[[1, end]]),
     wL2_12 = 25 * 12^2 / 12)
end

begin
    f = fig_member(hm, haunched; title = "haunched girder (deep ends): N, V, M")
    savefig("A5_haunched.png", f)
end

#=
## Numbers for the report

`reports/A5_joints_springs_variable.pdf` prints these (see `figures.jl`, `report_data`).
=#
begin
    # exact midspan moments from the member's force state (the sampled stations miss x = 6 m)
    fh = InternalForces(hm, haunched; resolution = 60)
    useMz = maximum(abs, fh.Mz) > maximum(abs, fh.My)
    sh = internal_forces(hm, haunched); seg, τ = Asap.locate_segment(haunched, 0.5)
    Mh_mid = useMz ? moment_z(sh[seg], τ) : moment_y(sh[seg], τ)
    sp = internal_forces(pm, prismatic)
    Mp_mid = useMz ? moment_z(sp, 0.5) : moment_y(sp, 0.5)
    report_data("A5";
        E = steel.E,
        column = (A = column.A, Ix = column.Ix), girder = (A = girder.A, Ix = girder.Ix),
        deep = (A = deep.A, Ix = deep.Ix),
        joints = (span = 8.0, height = 4.0, w = 20.0, H = 15.0, ks = ks, drift_mm = drifts, rigid_mm = rigid,
                  gap_k1e8 = (drifts[end] - rigid) / rigid),
        footing = (h = 4.0, H = 10.0, V = 50.0, k_trans = 1e6, EI = EI, ktheta = kθs, asap_mm = tops, exact_mm = exact.(kθs),
                   max_rel_err = maximum(abs.(tops .- exact.(kθs)) ./ exact.(kθs)),
                   flexural_mm = 1000 * 10 * 4^3 / (3EI),
                   translational_mm = 1000 * 10 / 1e6,                       # 10 kN on the 1e6 kN/m soil spring
                   max_abs_diff_mm = maximum(abs.(tops .- exact.(kθs)))),
        haunch = (span = 12.0, w = 25.0, fractions = [0.2, 0.8], wL2_8 = 25 * 12^2 / 8,
                  haunched_support = maximum(abs, Mh[[1, end]]), prismatic_support = maximum(abs, Mp[[1, end]]),
                  haunched_mid = Mh_mid, prismatic_mid = Mp_mid,
                  wL2_12 = 25 * 12^2 / 12, wL2_24 = 25 * 12^2 / 24))
end
