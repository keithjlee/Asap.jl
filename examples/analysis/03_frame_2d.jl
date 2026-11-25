#=
Walkthrough A3 — a 2D portal frame: loads, releases, internal forces

Shift+enter through the blocks. Units: kN, m.

A single-bay portal (8 m span, 4 m columns) in the XZ plane: a distributed
roof load, a point load, a lateral load, member force diagrams, the exact
deflected shape, and what end releases do to the moments.
=#

begin
    using Asap, LinearAlgebra
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

#=
## Sections

`Section(material, A, Ix, Iy, J)` — strong-axis `Ix`, weak-axis `Iy`.
=#
begin
    steel = Material(200e6, 77e6, 7.85, 0.3)
    column = Section(steel, 1.3e-2, 3.0e-4, 1.0e-4, 1.0e-6)     # ~W310×97
    girder = Section(steel, 1.1e-2, 4.5e-4, 1.2e-4, 1.0e-6)     # ~W460×82
end

#=
## The frame

Fixed bases; `planarize!(model, :ZX)` keeps it in its plane. Frame element
orientation: the roll angle sets which section axis bends in which plane —
with the default (π/2) these columns and girders bend about their STRONG
axis (`Ix`) in the XZ plane; `rollangle = 0` would turn them to the weak axis.
=#
begin
    function portal(; release = :fixedfixed)
        a = Node([0.0, 0.0, 0.0], :fixed); b = Node([0.0, 0.0, 4.0], :free)
        c = Node([8.0, 0.0, 4.0], :free); d = Node([8.0, 0.0, 0.0], :fixed)
        colL = FrameElement(a, b, column, :column)
        beam = FrameElement(b, c, girder, :girder; release = release)
        colR = FrameElement(d, c, column, :column)
        loads = AbstractLoad{Float64}[LineLoad(beam, [0.0, 0.0, -20.0]),          # 20 kN/m roof load
                                      PointLoad(beam, 0.25, [0.0, 0.0, -30.0]),    # 30 kN at quarter span
                                      NodeForce(b, [15.0, 0.0, 0.0])]              # wind at eave
        m = Model([a, b, c, d], AbstractElement{Float64}[colL, beam, colR], loads)
        planarize!(m, :ZX)
        solve!(m)
        return m, beam, colL
    end
    frame, beam, colL = portal()
    res = frame.results
    (drift_mm = 1000 * displacement(res, beam.nodeStart)[1], base_moment = reaction(res, colL.nodeStart)[5])
end

#=
## Internal forces along a member

`InternalForces(model, element)` samples N, V, M along the member (exact,
including the jump at the point load). Peak girder moment:
=#
begin
    fi = InternalForces(frame, beam; resolution = 50)
    # the bending component that is active for this orientation
    M = maximum(abs, fi.Mz) > maximum(abs, fi.My) ? fi.Mz : fi.My
    (max_moment = maximum(M), min_moment = minimum(M))
end

begin
    f = fig_member(frame, beam; title = "girder", width = 210)
    savefig("A3_girder_diagrams.pdf", f)
end

begin
    diag = maximum(abs, InternalForces(frame, beam).Mz) > maximum(abs, InternalForces(frame, beam).My) ? :Mz : :My
    f = fig_lines(frame; plane = :XZ, scale = 50.0, diagram = diag, color = kjl_darkgray,
        title = "moment, deflection ×50", width = W_HALF)
    savefig("A3_portal.pdf", f)
end

#=
## End releases

`release = :freefree` pins the girder at both ends: it becomes simply
supported (midspan moment → wL²/8 = 160 kN·m from the UDL alone) and the
columns become cantilevers, so the frame drifts far more.
=#
begin
    pinned, beam_p, _ = portal(; release = :freefree)
    fp = InternalForces(pinned, beam_p)
    Mp = maximum(abs, fp.Mz) > maximum(abs, fp.My) ? fp.Mz : fp.My
    (fixed_drift_mm = 1000 * displacement(res, beam.nodeStart)[1],
     pinned_drift_mm = 1000 * displacement(pinned.results, beam_p.nodeStart)[1],
     pinned_girder_max_moment = maximum(abs, Mp))
end

#=
## Numbers for the report

`reports/A3_frame_2d.pdf` prints these (see `figures.jl`, `report_data`).
Hand checks: statics of the reactions; the pinned girder is a simply supported
beam (UDL w plus P at L/4), and its columns are two equal cantilevers sharing
the eave load.
=#
begin
    Lg, Hc, w, P, a_P, Hw = 8.0, 4.0, 20.0, 30.0, 0.25 * 8.0, 15.0
    RA = w * Lg / 2 + P * (1 - a_P / Lg)                          # simply supported left reaction
    xstar = max((RA - P) / w, a_P)                                # zero shear (right of the point load)
    ss_max = RA * xstar - w * xstar^2 / 2 - P * (xstar - a_P)
    supports = (colL.nodeStart, frame.elements[3].nodeStart)
    rxn(r, n) = (Fx = reaction(r, n)[1], Fz = reaction(r, n)[3], My = reaction(r, n)[5])
    psupports = (pinned.elements[1].nodeStart, pinned.elements[3].nodeStart)
    report_data("A3";
        span = Lg, height = Hc, w = w, P = P, a_P = a_P, H = Hw, E = steel.E,
        column = (A = column.A, Ix = column.Ix, Iy = column.Iy), girder = (A = girder.A, Ix = girder.Ix, Iy = girder.Iy),
        drift_mm = 1000 * displacement(res, beam.nodeStart)[1],
        reactions = [rxn(res, n) for n in supports],
        sum_Fx = sum(n -> reaction(res, n)[1], supports), sum_Fz = sum(n -> reaction(res, n)[3], supports),
        total_vertical = w * Lg + P,
        girder_max_moment = maximum(M), girder_min_moment = minimum(M),
        moment_component = maximum(abs, fi.Mz) > maximum(abs, fi.My) ? "Mz" : "My",
        pinned = (drift_mm = 1000 * displacement(pinned.results, beam_p.nodeStart)[1],
                  girder_max_moment = maximum(abs, Mp),
                  reactions = [rxn(pinned.results, n) for n in psupports],
                  RA_hand = RA, ss_max_hand = ss_max, girder_moment_at_xstar = let st = internal_forces(pinned, beam_p)
                      maximum(abs, fp.Mz) > maximum(abs, fp.My) ? moment_z(st, xstar / Lg) : moment_y(st, xstar / Lg)
                  end, ss_x_hand = xstar, wL2_8 = w * Lg^2 / 8,
                  drift_cantilever_hand_mm = 1000 * (Hw / 2) * Hc^3 / (3 * steel.E * column.Ix),
                  base_moment_hand = Hw / 2 * Hc))
end
