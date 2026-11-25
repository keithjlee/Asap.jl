#=
Walkthrough 7 — a dapped-end beam: cuts, ties and a singular corner

Shift+enter through the `begin … end` blocks top to bottom (or run the file).
Units: kN, m.

A precast beam (6 m long, 900 mm deep, 300 mm wide) is dapped at both
ends: the bottom half of the section is cut away over the last 300 mm, and
the beam sits on its nibs (450 mm deep) on bearing pads. The dap is the
textbook D-region: the reaction enters a short nib, which must hand it to
the full-depth beam around a re-entrant corner. A strut-and-tie model of it
has a horizontal tie along the bottom of the nib and a vertical hanger tie
just past the dap face, and cracks start at the re-entrant corner.

What this shows, beyond walkthrough 6:
- interior `segments` in `mesh_polygon`: element edges along chosen lines,
  so `section_cut` can integrate the exact force across them;
- reactions applied as pad pressures (a flexible bearing), with the beam
  held only against rigid-body motion: a self-equilibrated load set;
- the tie forces of a strut-and-tie model read off as integrated cut forces,
  which converge under refinement, while the stress at the re-entrant
  corner does not (it is singular).

Environment: `examples/Project.toml` (Asap and AsapOptim by path + plotting +
DelaunayTriangulation).
    julia --project=examples examples/shells/07_dapped_end.jl
=#

begin
    using Asap, LinearAlgebra
    using DelaunayTriangulation                        # activates mesh_polygon
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

#=
## Geometry

Plane coordinates (x along the beam, z up), placed in the global XZ plane.
The outline runs counter-clockwise; every point where a pad or the floor load
starts or ends is a vertex. The floor load acts on the full-depth part only
(from 0.6 m to 5.4 m), so no load sits on the nib or on the cut lines.

Interior `segments` are lines the element edges must follow:
- the nib/beam interface, the vertical through the re-entrant corner (x = 0.3);
- the dap-top level (z = 0.45) from the corner into the span, where the
  hanger tension crosses;
- a full-depth section one beam depth past the dap face (x = 1.2).
Segments may meet (split them where they cross) and may end on the outline.
=#
begin
    L, H, hd, ln, t = 6.0, 0.9, 0.45, 0.30, 0.30       # length, depth, dap height, nib length, width
    pad = (0.05, 0.20)                                 # pad under each nib, 150 mm
    xb = sum(pad) / 2                                  # pad centre
    brgL = ((pad[1], hd), (pad[2], hd))
    brgR = ((L - pad[2], hd), (L - pad[1], hd))
    floor_load = ((2ln, H), (L - 2ln, H))
    outline = [(ln, 0.0), (L - ln, 0.0), (L - ln, hd), brgR..., (L, hd), (L, H), reverse(floor_load)...,
               (0.0, H), (0.0, hd), brgL..., (ln, hd)]
    xc = ln + H                                        # one depth past the dap face
    interface = ((ln, hd), (ln, H))
    daptop = (((ln, hd), (xc, hd)), ((xc, hd), (L / 2, hd)))
    section = (((xc, 0.0), (xc, hd)), ((xc, hd), (xc, H)))
    concrete = Material(30.0e6, 2.4, 0.2)               # E [kN/m²], ρ [t/m³], ν
    beam_sec = ShellSection(concrete, t)
end

function beam_mesh(h)
    mesh_polygon(outline, beam_sec; segments = [interface, daptop..., section...],
        points = [(xb, hd), (L - xb, hd)],             # pad centres: the rigid-body supports
        h = h, min_angle = 28, e1 = [1.0, 0, 0], e2 = [0.0, 0, 1], id = :beam)
end

begin
    h = 0.05
    mesh = beam_mesh(h)
    (nodes = length(mesh.nodes), triangles = length(mesh.elements),
     interface_nodes = length(nodes_on_segment(mesh, interface...)))
end

#=
## Loads and supports: a self-equilibrated set

A bearing pad is flexible: it presses on the nib with a roughly uniform
pressure, and it cannot pull. Holding the pad's nodes fixed would instead
make the pad rigid and bilateral: the nib would rotate on it, the pressure
would pile up at the pad's inner edge and turn into tension at its outer
edge. So the reactions enter as LOADS: a uniform upward pressure over each
pad, equal to the reaction the statics of the whole beam give (R = 300 kN
per end, the floor load w = 2R / 4.8 m = 125 kN/m).

The loads then balance by themselves. The beam is held only against
rigid-body motion: one node at the left pad centre pinned, one at the right
pad centre on a roller (out-of-plane DOFs held everywhere, as in walkthrough
6). Those supports carry nothing; their reactions measure how exactly the
load set balances.
=#
function beam_model(mesh)
    for n in mesh.nodes
        fixnode!(n, Bool[1, 0, 1, 0, 1, 0])           # in-plane panel (XZ)
    end
    at(x, z) = only(filter(n -> hypot(n.position[1] - x, n.position[3] - z) < 1e-9, mesh.nodes))
    pin, roller = at(xb, hd), at(L - xb, hd)
    fixnode!(pin, Bool[1, 1, 1, 1, 1, 0])
    fixnode!(roller, Bool[0, 1, 1, 1, 1, 0])
    R = 300.0
    w = 2R / (floor_load[2][1] - floor_load[1][1])
    p = R / (pad[2] - pad[1])
    loads = AbstractLoad{Float64}[edge_loads(mesh, floor_load..., [0.0, 0, -w])...,
        edge_loads(mesh, brgL..., [0.0, 0, p])..., edge_loads(mesh, brgR..., [0.0, 0, p])...]
    model = Model(mesh.nodes, AbstractElement{Float64}[mesh.elements...], loads)
    solve!(model)
    return model, R, w, pin, roller
end

begin
    model, R, w, pin, roller = beam_model(mesh)
    support = maximum(abs, [reaction(model.results, pin); reaction(model.results, roller)])
    (R = R, w = w, support_reactions = support)
end

#=
## Cut forces

`section_cut(model, elements, cut_nodes)` sums the elements' nodal forces at
the nodes of a cut: the force and moment one side transmits to the other,
EXACT for equilibrium on any mesh (no stress smoothing). Two kinds of use:

- the whole cut (all its nodes): the section's resultants V, N, M;
- the cut's nodes one at a time, accumulated along it: the force profile
  across the cut, whose tension part is a strut-and-tie model's tie.

The nib/beam interface (the free body is the nib, x < 0.3, above the dap)
must carry V = R and M = R·(0.3 − 0.125) = 52.5 kN·m about the re-entrant
corner, whatever the mesh.
=#
# force through the nodes of a cut, accumulated along it (nodes in order), in direction `dir`
cut_profile(model, els, nodes, dir) = cumsum([dot(section_cut(model, els, [n])[1], dir) for n in nodes])

function cuts(model, mesh)
    nib = elements_in(mesh, [(0.0, hd), (ln, hd), (ln, H), (0.0, H)])
    inode = nodes_on_segment(mesh, interface...)                  # corner → top
    F, M = section_cut(model, nib, inode; about = [ln, 0.0, hd])
    # nib tie: horizontal tension across the interface, from the corner up
    T = cut_profile(model, nib, inode, [1.0, 0, 0])
    kt = argmax(T)
    # hanger: vertical tension across the dap-top level, from the corner into the
    # span. The free body below it (the bottom strip) carries no load, so the total
    # over the whole level is zero: tension next to the corner, compression further in
    strip = elements_in(mesh, [(ln, 0.0), (L - ln, 0.0), (L - ln, hd), (ln, hd)])
    hnode = [nodes_on_segment(mesh, daptop[1]...); nodes_on_segment(mesh, daptop[2]...)[2:end]]
    Hz = cut_profile(model, strip, hnode, [0.0, 0, 1.0])
    kh = argmax(Hz)
    # the full-depth section at x = 1.2 (free body: everything left of it)
    left = elements_in(mesh, [(0.0, hd), (ln, hd), (ln, 0.0), (xc, 0.0), (xc, H), (0.0, H)])
    snode = [nodes_on_segment(mesh, section[1]...); nodes_on_segment(mesh, section[2]...)[2:end]]
    _, Ms = section_cut(model, left, snode; about = [xc, 0.0, H / 2])
    Sx = cut_profile(model, left, snode, [1.0, 0, 0])
    # the three profiles (position along the cut, accumulated force), for the figure
    profiles = (nib = (s = [n.position[3] - hd for n in inode], F = T),
                section = (s = [n.position[3] for n in snode], F = Sx),
                hanger = (s = [n.position[1] - ln for n in hnode], F = Hz))
    return (V = -F[3], M = -M[2], nib_tie = T[kt], nib_tie_depth = inode[kt].position[3] - hd,
            hanger = Hz[kh], hanger_width = hnode[kh].position[1] - ln,
            Ms = -Ms[2], section_tension = maximum(Sx), profiles)
end

begin
    c = cuts(model, mesh)
    (V = c.V, M = c.M, M_statics = R * (ln - xb))
end

#=
## The nib tie

The interface's moment is carried by a tension in the lower part of the nib
and a compression above it. The tension, integrated from the corner up to
where it turns into compression, is the force the nib's horizontal
reinforcement takes in a strut-and-tie model. Its lever arm z = M / T is the
internal lever of the nib.
=#
begin
    (nib_tie = round(c.nib_tie, digits = 1), tension_zone_depth = round(c.nib_tie_depth, digits = 3),
     lever = round(c.M / c.nib_tie, digits = 3), lever_over_nib_depth = round(c.M / c.nib_tie / (H - hd), digits = 2))
end

#=
## The hanger

The reaction has to get from the nib into the full-depth beam. A
strut-and-tie model hangs it up with a vertical tie just past the dap face:
the hanger, sized for the full reaction R. The elastic, UNCRACKED field
shares the transfer differently: the vertical tension crossing the dap-top
level next to the corner integrates to about two thirds of R, over a width
of about a quarter of a metre. Do not size a hanger from that figure: the
corner cracks at a fraction of this load (the elastic stress there, N1/t, is
several times the concrete's tensile strength, and grows with refinement),
and after that the load path is the cracked one the strut-and-tie model
describes: the hanger needs the full R. The elastic field tells where the hanger
goes (next to the dap face) and how wide its zone is.
=#
begin
    (hanger = round(c.hanger, digits = 1), hanger_over_R = round(c.hanger / R, digits = 2),
     hanger_zone_width = round(c.hanger_width, digits = 3))
end

#=
## One depth past the dap: a beam again

At x = 1.2 m (one beam depth past the dap face) the section is a B-region: the
cut's moment is M = R·(1.2 − 0.125) − w·0.6²/2 by statics, and its tension
resultant matches that of the linear (beam-theory) stress block, M / (2H/3).
(The cut's shear differs from R − w·0.6 by the share of the floor load sitting
on the cut's top node, w·s/2 for an element edge s there: a load applied AT a
cut node belongs to neither side. The interface cut above avoids that, as no
load touches it.)
=#
begin
    Ms_statics = R * (xc - xb) - w * (xc - floor_load[1][1])^2 / 2
    T_beam = Ms_statics / (2H / 3)
    (Ms = round(c.Ms, digits = 3), Ms_statics = Ms_statics,
     section_tension = round(c.section_tension, digits = 1), beam_theory = round(T_beam, digits = 1))
end

#=
## The corner: where cracks start

Principal tension peaks at the re-entrant corner. Its direction there is
inclined at about 45° (averaged over a ring 50–120 mm from the corner, inside
the full-depth part), so the crack that starts at the corner runs at about
45° up into the beam, perpendicular to it. That crack is why the dap needs
the hanger (across it) and why codes ask for inclined or well-anchored bars
at the corner.
=#
begin
    g = ShellGeo(model)
    ring = [k for k in axes(g.centers, 2)
            if 0.05 < hypot(g.centers[1, k] - ln, g.centers[3, k] - hd) < 0.12 &&
               g.centers[1, k] > ln && g.centers[3, k] > hd]
    # N1 direction in global axes (θ is measured from each element's own e₁), averaged
    # as a direction field (doubled angle), weighted by N1
    c2 = sum(ring) do k
        Rf, θ = g.frames[k], g.θ[k]
        d1 = cos(θ) * Rf[1, :] + sin(θ) * Rf[2, :]
        g.N1[k] * cis(2atan(d1[3], d1[1]))
    end
    tension_dir = mod(rad2deg(angle(c2) / 2), 180)      # from +x, counter-clockwise
    crack_dir = mod(tension_dir - 90, 180)
    (tension_direction = round(tension_dir, digits = 1), crack_direction = round(crack_dir, digits = 1),
     peak_N1 = round(maximum(g.N1), digits = 0))
end

#=
## Figures

The whole beam, then a close-up of the left end (the dap is too small to read at full length)
with the three cut lines, and the force profiles along the cuts: the peak of each is a tie.
=#
begin
    win = ((0.0, 1.5), (0.0, H))                       # the left end, past the full-depth section
    cutlines = [(seg = interface, label = "interface"),
                (seg = ((ln, hd), (L / 2, hd)), label = "dap top", at = (0.95, hd),
                 align = (:center, :bottom), offset = (0, 3)),
                (seg = ((xc, H), (xc, 0.0)), label = "x = 1.2")]
    savefig("S7_beam.pdf", fig_principal(model, mesh; geo = g, title = "principal membrane forces",
        bearings = [brgL, brgR], loads = [floor_load], width = W_FULL))
    savefig("S7_principal.pdf", fig_principal(model, mesh; geo = g, title = "left end",
        bearings = [brgL], loads = [floor_load], window = win, cuts = cutlines, width = W_WIDE))
    savefig("S7_N1.pdf", fig_wall(model, mesh; geo = g, field = :N1, title = "N_1 tension, left end",
        window = win, cuts = cutlines, width = W_WIDE))
    p = c.profiles
    savefig("S7_cuts.pdf", fig_cutprofiles([
        (s = p.nib.s, F = p.nib.F, vertical = true, title = "nib interface",
         slabel = "z above dap [m]", flabel = "F_x [kN]"),
        (s = p.section.s, F = p.section.F, vertical = true, title = "section x = 1.2",
         slabel = "z [m]", flabel = "F_x [kN]"),
        (s = p.hanger.s, F = p.hanger.F, title = "dap top",
         slabel = "x past dap face [m]", flabel = "F_z [kN]")]))
end

#=
## Mesh study

The cut forces converge as the mesh is refined (the interface's V and M
are exact on every mesh); the peak tension at the corner keeps growing: the
elastic field is singular at a re-entrant corner. A strut-and-tie model reads
integrated forces and directions, never the corner peak.
=#
begin
    study = map((0.1, 0.05, 0.025)) do hh
        m = beam_mesh(hh)
        md, = beam_model(m)
        ch = cuts(md, m)
        (h = hh, triangles = length(m.elements), nib_tie = ch.nib_tie, hanger = ch.hanger,
         peak_N1 = maximum(ShellGeo(md).N1))
    end
    study
end

#=
## Numbers for the report

`reports/S7_dapped_end.pdf` prints these (see `figures.jl`, `report_data`).
=#
report_data("S7";
    L, H, hd, ln, t, pad_start = pad[1], pad_end = pad[2], xb, R, w, h,
    nodes = length(mesh.nodes), triangles = length(mesh.elements), support_reactions = support,
    V = c.V, M = c.M, M_statics = R * (ln - xb), nib_tie = c.nib_tie, nib_tie_depth = c.nib_tie_depth,
    nib_lever = c.M / c.nib_tie, hanger = c.hanger, hanger_width = c.hanger_width,
    Ms = c.Ms, Ms_statics, section_tension = c.section_tension, T_beam,
    tension_dir, crack_dir, peak_N1 = maximum(g.N1),
    study_h = [s.h for s in study], study_triangles = [s.triangles for s in study],
    study_nib_tie = [s.nib_tie for s in study], study_hanger = [s.hanger for s in study],
    study_peak = [s.peak_N1 for s in study])
