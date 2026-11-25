#=
Walkthrough 6 — a disturbed region: a deep transfer panel with openings

Shift+enter through the `begin … end` blocks top to bottom (or run the file).
Units: kN, m.

A transfer panel (8 m long, 4 m deep, 400 mm thick) carries a column load
on its top edge down to two bearing pads, around a circular duct and a
rectangular opening. Nothing about it is beam-like: it is a D-region, the
regime a strut-and-tie model is built for. The linear-elastic stress field
of the panel is the usual guide for that model — where the compression
struts run, where the ties must be, how the load splits between the
bearings.

What this shows: meshing an arbitrary polygon with holes into triangular
shells (`mesh_polygon`, through DelaunayTriangulation.jl), placing bearings
and edge loads on boundary segments, reading reactions per bearing, and the
principal membrane forces that sketch the strut-and-tie model. A mesh study
closes it: reactions converge immediately, peak stresses at re-entrant
corners do not (they are singular), which is why a strut-and-tie model reads
the field's directions, not its peaks.

Environment: `examples/Project.toml` (Asap and AsapOptim by path + plotting +
DelaunayTriangulation).
    julia --project=examples examples/shells/06_d_region.jl
=#

begin
    using Asap, LinearAlgebra, Statistics
    using DelaunayTriangulation                        # activates mesh_polygon
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

#=
## Geometry

Plane coordinates (x along the panel, z up), placed in the global XZ plane
(`e1 = X`, `e2 = Z`). Every point where a bearing or the column load starts
or ends is a vertex of the outline, so the mesh has nodes exactly there.
=#
begin
    L, H, t = 8.0, 4.0, 0.40
    brgL = ((0.4, 0.0), (1.2, 0.0))                    # bearing pads, 0.8 m each
    brgR = ((6.8, 0.0), (7.6, 0.0))
    col = ((3.6, H), (4.4, H))                         # column above, 0.8 m wide
    outline = [(0.0, 0.0), brgL..., brgR..., (L, 0.0), (L, H), reverse(col)..., (0.0, H)]
    duct = circle_polygon((2.3, 2.2), 0.45; n = 36)     # Ø 900 mm duct
    opening = [(5.3, 0.9), (6.3, 0.9), (6.3, 2.5), (5.3, 2.5)]
    concrete = Material(30.0e6, 2.4, 0.2)               # E [kN/m²], ρ [t/m³], ν
    panel_sec = ShellSection(concrete, t)
end

#=
## Mesh

`mesh_polygon` triangulates the outline minus the holes (constrained
Delaunay, refined to triangles of edge ≈ `h` with every angle ≥ `min_angle`).
The elements are ANDES-OPT/MITC3+ triangles; the mesh keeps the nodes of
every input edge (`boundary_nodes`) and answers segment queries.
=#
function panel_mesh(h)
    mesh_polygon(outline, panel_sec; holes = [duct, opening], h = h, min_angle = 28,
        e1 = [1.0, 0, 0], e2 = [0.0, 0, 1], id = :panel)
end

begin
    h = 0.15
    mesh = panel_mesh(h)
    min_angle = minimum(mesh.elements) do el
        x = map(n -> n.position, el.nodes)
        minimum(k -> acosd(clamp(dot(normalize(x[mod1(k + 1, 3)] - x[k]), normalize(x[mod1(k + 2, 3)] - x[k])), -1, 1)), 1:3)
    end
    (nodes = length(mesh.nodes), triangles = length(mesh.elements), min_angle = round(min_angle, digits = 1),
     area = round(sum(Asap.area, mesh.elements), digits = 4))
end

#=
## Supports and loads

The panel works in its own plane (XZ): every node is held out of plane (Y)
and against rotation about the in-plane axes (X, Z); the drilling rotation
(about Y) stays free. Along each bearing pad the translations are held, and —
because the material along a held line cannot rotate — the drilling rotation
too. The left pad is pinned, the right one a roller (free in X), so the
panel is statically determinate as a body.

The column delivers 2400 kN as a uniform line load over its 0.8 m width
(`edge_loads` on that boundary segment adds the drilling moments a drilling
membrane needs for a uniform traction). Self-weight acts as well.
=#
function panel_model(mesh)
    for n in mesh.nodes
        fixnode!(n, Bool[1, 0, 1, 0, 1, 0])
    end
    foreach(n -> fixnode!(n, Bool[0, 0, 0, 0, 0, 0]), nodes_on_segment(mesh, brgL...))
    foreach(n -> fixnode!(n, Bool[1, 0, 0, 0, 0, 0]), nodes_on_segment(mesh, brgR...))
    P = 2400.0
    w = P / 0.8
    loads = AbstractLoad{Float64}[edge_loads(mesh, col..., [0.0, 0.0, -w])...,
        [SelfWeight(el; g = [0.0, 0, -9.81]) for el in mesh.elements]...]
    model = Model(mesh.nodes, AbstractElement{Float64}[mesh.elements...], loads)
    solve!(model)
    return model, P
end

begin
    model, P = panel_model(mesh)
    W = areal_mass(panel_sec) * 9.81 * sum(Asap.area, mesh.elements)
    nL, nR = nodes_on_segment(mesh, brgL...), nodes_on_segment(mesh, brgR...)
    FL, ML = total_reaction(model.results, nL; about = zeros(3))
    FR, MR = total_reaction(model.results, nR; about = zeros(3))
    (left = round(FL[3], digits = 1), right = round(FR[3], digits = 1), total = round(FL[3] + FR[3], digits = 1),
     applied = round(P + W, digits = 1))
end

#=
## Statics of the reactions

The reactions must balance the loads exactly (vertical force and moment
about the origin), whatever the mesh. How the load SPLITS between the pads
depends on where each pad's resultant acts: a lever rule that puts it at the
pad centre is only approximate, because the bearing pressure is not uniform —
it peaks at the pad's inner edge. The resultant's actual position, from the
reaction moment, is where a strut-and-tie model puts its bottom node. (The
self-weight's centroid is off midspan too: the openings remove material
unevenly.)
=#
begin
    xc = sum(el -> Asap.area(el) * sum(n -> n.position[1], el.nodes) / 3, mesh.elements) / sum(Asap.area, mesh.elements)
    xL, xR = (brgL[1][1] + brgL[2][1]) / 2, (brgR[1][1] + brgR[2][1]) / 2
    # the reaction resultant of each pad acts at its centre only approximately (the
    # distribution under a pad is not uniform), so use the exact moment balance:
    moment_residual = ML[2] + MR[2] + (4.0 * P + xc * W)     # ≈ 0 about the Y axis
    R_right_lever = (P * (4.0 - xL) + W * (xc - xL)) / (xR - xL)
    # where each pad's resultant acts (pads at z = 0: M_y = −x·F_z, couples included)
    xL_eff, xR_eff = -ML[2] / FL[3], -MR[2] / FR[3]
    (moment_residual = moment_residual, right_by_pad_centres = round(R_right_lever, digits = 1),
     right_fe = round(FR[3], digits = 1), left_resultant_at = round(xL_eff, digits = 3),
     right_resultant_at = round(xR_eff, digits = 3))
end

#=
## The stress field: struts and ties

Principal membrane forces per element (`ShellGeo`: N1 ≥ N2, angle θ from the
element's own local axis). Compression (N2 < 0) runs from the column down to the
bearings, around the openings; tension (N1 > 0) concentrates along the
bottom between the pads — the tie. Between each pad and the column the
compression fans out from the column and gathers into the pad: a
bottle-shaped strut. Its local direction part-way up is steeper than the
straight line from the pad to the column, which a strut-and-tie model draws
as the strut's resultant; the spreading is what the model's transverse ties
(or the bottle-strut rules of the design codes) account for.
=#
begin
    g = ShellGeo(model)
    near(x, z) = argmin(k -> hypot(g.centers[1, k] - x, g.centers[3, k] - z), axes(g.centers, 2))
    # strut angle above the left pad: principal compression direction vs the line to the column
    # The strut from the left pad to the column: sample a quarter of the way up the
    # line pad centre (0.8, 0) → column centre (4, 4), clear of the duct. θ is measured
    # from each element's OWN e₁ axis (for triangles: its first edge), so rebuild the N2
    # direction in global axes from the element frame, and average it over a small
    # disk weighted by |N2| (a direction field: average the doubled angle)
    sx, sz = 0.8 + 0.25 * 3.2, 0.25 * 4.0
    disk = [k for k in axes(g.centers, 2) if hypot(g.centers[1, k] - sx, g.centers[3, k] - sz) < 0.3]
    c2 = sum(disk) do k
        R, θ = g.frames[k], g.θ[k]
        d2 = -sin(θ) * R[1, :] + cos(θ) * R[2, :]        # N2 ⟂ N1, in the panel plane
        φ = atan(d2[3], d2[1])
        abs(g.N2[k]) * cis(2φ)
    end
    strut_fe = mod(rad2deg(angle(c2) / 2), 180)
    strut_geom = rad2deg(atan(4.0, 3.2))
    # tie: tension along the bottom at midspan, and its height
    bottom = [k for k in axes(g.centers, 2) if abs(g.centers[1, k] - 4.0) < 0.3]
    kt = argmax(k -> g.N1[k], bottom)
    (strut_angle_fe = round(strut_fe, digits = 1), strut_angle_load_to_pad = round(strut_geom, digits = 1),
     peak_tension_midspan = round(g.N1[kt], digits = 1), at_z = round(g.centers[3, kt], digits = 2),
     peak_compression = round(minimum(g.N2), digits = 1))
end

begin
    savefig("S6_principal.pdf", fig_principal(model, mesh; geo = g, title = "principal membrane forces",
        bearings = [brgL, brgR], loads = [col]))
    savefig("S6_N2.pdf", fig_wall(model, mesh; geo = g, field = :N2, title = "N_2 compression [kN/m]"))
end

#=
## Mesh study

Reactions are integrals of the field and converge at once. Peak principal
forces sit at re-entrant corners (the rectangular opening, the pad edges),
where the elastic field is singular: they grow without bound under
refinement. Read a strut-and-tie model from the field's DIRECTIONS and its
resultants, never from corner peaks.
=#
begin
    study = map((0.4, 0.2, 0.1)) do h
        m = panel_mesh(h)
        mod_, _ = panel_model(m)
        FRh, _ = total_reaction(mod_.results, nodes_on_segment(m, brgR...))
        gh = ShellGeo(mod_)
        (h = h, triangles = length(m.elements), right_reaction = FRh[3], peak_compression = minimum(gh.N2))
    end
    study
end

#=
## Numbers for the report

`reports/S6_d_region.pdf` prints these (see `figures.jl`, `report_data`).
=#
report_data("S6";
    L, H, t, P, W, h, nodes = length(mesh.nodes), triangles = length(mesh.elements), min_angle,
    left = FL[3], right = FR[3], moment_residual, right_by_pad_centres = R_right_lever,
    left_resultant_at = xL_eff, right_resultant_at = xR_eff,
    strut_fe, strut_geom, peak_tension = g.N1[kt], peak_tension_z = g.centers[3, kt],
    peak_compression = minimum(g.N2),
    study_h = [s.h for s in study], study_triangles = [s.triangles for s in study],
    study_right = [s.right_reaction for s in study], study_peak = [s.peak_compression for s in study])
