#=
Figure generators for the walkthroughs.

Plain functions of Asap's plot-extraction objects (`ShellGeo`, `WallMesh`, the solved
`Model`) that return Makie figures in the kjlMakie house style:

    fig_wall(model, w; field = :N2, scale = 0.0)    wall elevation, contoured field on the (scaled) deformed mesh
    fig_principal(model, w)                         principal membrane-force trajectories (ties / struts)
    fig_bars(labels, values; ylabel)                bar chart (pier shears, storey shears, …)
    fig_slab(model; field = :Mxx, elements)         plan view of a slab, contoured
    fig_3d(model; field = :N2, scale)               3D view of shells + frames, deformed
    fig_series(xs, series; xlabel, ylabel, ref)     convergence / parametric curves
    fig_lines(model; plane, scale, diagram)         2D view of frames/trusses: deformed (exact member curves),
                                                    coloured by axial force, optional moment diagram
    fig_member(model, el)                           N / V / M diagrams along one member

House style (kjlMakie; the rules are in AsapEcosystem/CLAUDE.md "Figures"):
- internal-force diagrams: ONE diagram black with a 50 % fill (`diagramband!`); several compared
  `kjl_categorical` (blue, green, yellow, grays).
- colours by ROLE only: kjl_tension (blue, +) / kjl_compression (red, −), signed fields
  kjl_diverging (red → white → blue, `symmetric_range`), neutral fields such as displacement
  `neutral_colormap` (white → orange), one highlighted series kjl_emphasis (blue) with the rest
  grayscale (`kjl_series(i)`), a bad item kjl_alert (red). No literal hues, no Makie colormaps.
- data axes: theme `kjl_light`; structure drawings: `structurestyle!` + `fitstructure!`
  (true aspect, padded limits, no axes).
- labels in CAPS (`savefig` applies `capslabels!`; units in [brackets], Greek, rich-text
  subscripts and `{…}` keep their case).
- figure units are points: every generator takes `width` = the PRINTED width [pt]
  (W_FULL / W_WIDE / W_PLAN / W_HALF below); `savefig` writes a PDF at 1 pt per unit and the
  reports embed it at natural size, so the 9 pt text prints at 9 pt.

`include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)` from any
walkthrough (it loads the plotting stack); figures go to `figures/` next to the walkthrough.
`report_data(name; key = value, …)` writes the numbers a walkthrough's report prints
to `reports-src/data/<name>.json` (the report sources are not distributed; the PDFs are in
`reports/`).
=#

using Asap, LinearAlgebra
using kjlMakie, CairoMakie
CairoMakie.activate!()
set_theme!(kjl_light)

# printed widths [pt] (kjl-typst report: 504 pt text width)
const W_FULL = 470                # full width (a little inside the text block)
const W_WIDE = 428                # 85 %
const W_PLAN = 353                # 70 %
const W_HALF = 240                # one column of a two-column grid

# each walkthrough saves into `figures/` next to itself: call `set_figdir(@__DIR__)`
const FIGDIR = Ref(joinpath(@__DIR__, "figures"))
set_figdir(dir) = (FIGDIR[] = joinpath(dir, "figures"); mkpath(FIGDIR[]); FIGDIR[])

"""
    savefig(name, fig)

Caps the labels and writes `figures/<name>.pdf` (vector, 1 pt per figure unit: embed at
natural size). Any extension given in `name` is replaced by `.pdf`.
"""
function savefig(name, fig)
    capslabels!(fig)
    path = joinpath(FIGDIR[], splitext(name)[1] * ".pdf")
    savereport(path, fig)
    return fig
end

# ── report data ──────────────────────────────────────────────────────────────
# Each walkthrough ends with `report_data("A1"; key = value, …)`: the numbers its
# report (`reports/<id>_*.pdf`) prints, written to `reports-src/data/<name>.json`,
# so no number in a report is typed by hand. Values: numbers, strings, symbols,
# vectors and (named) tuples of those; non-finite floats are written as null.
using JSON
const REPORT_DATA = joinpath(@__DIR__, "reports-src", "data")
_jsonable(x::AbstractFloat) = isfinite(x) ? x : nothing
_jsonable(x::Symbol) = string(x)
_jsonable(x::Union{AbstractVector,Tuple}) = map(_jsonable, collect(x))
_jsonable(x::AbstractRange) = _jsonable(collect(x))
_jsonable(x::NamedTuple) = Dict(string(k) => _jsonable(v) for (k, v) in pairs(x))
_jsonable(x::AbstractDict) = Dict(string(k) => _jsonable(v) for (k, v) in x)
_jsonable(x) = x
function report_data(name::AbstractString; kw...)
    mkpath(REPORT_DATA)
    path = joinpath(REPORT_DATA, name * ".json")
    write(path, JSON.json(_jsonable(NamedTuple(kw)); pretty = true))
    return path
end

# ── fields ───────────────────────────────────────────────────────────────────

# colourbar labels (rich text keeps subscripts lowercase; units kN, m — displacements in mm)
_sub(base, s, rest) = rich(base, subscript(s), rest)
const FIELD_LABELS = Dict(
    :Nxx => _sub("N", "xx", "  [kN/m]"), :Nyy => _sub("N", "yy", "  [kN/m]"), :Nxy => _sub("N", "xy", "  [kN/m]"),
    :N1 => _sub("N", "1", "  MAJOR PRINCIPAL [kN/m]"), :N2 => _sub("N", "2", "  MINOR PRINCIPAL [kN/m]"),
    :Mxx => _sub("M", "xx", "  [kN·m/m]"), :Myy => _sub("M", "yy", "  [kN·m/m]"),
    :Mxy => _sub("M", "xy", "  [kN·m/m]"), :Qx => _sub("Q", "x", "  [kN/m]"), :Qy => _sub("Q", "y", "  [kN/m]"),
    :w => rich("w  OUT-OF-PLANE [mm]"), :umag => rich("|u|  [mm]"))

const DISPLACEMENT_FIELDS = (:w, :umag)
# displacement fields are drawn in mm
field_scale(field) = field in DISPLACEMENT_FIELDS ? 1000.0 : 1.0

"""
colormap and colorrange for `field` over `vals`: forces/moments are signed (kjl_diverging,
zero white at the centre); displacements are neutral (white → kjl_orange, zero white)
"""
field_colormap(field, vals) = field in DISPLACEMENT_FIELDS ? neutral_colormap(vals) : (kjl_diverging, symmetric_range(vals))

"per-shell-element values of `field` (`:w`, `:umag` are nodal and handled by `nodal_field`)"
function element_field(g::ShellGeo, field::Symbol)
    field === :Nxx && return g.N[1, :]
    field === :Nyy && return g.N[2, :]
    field === :Nxy && return g.N[3, :]
    field === :N1 && return g.N1
    field === :N2 && return g.N2
    field === :Mxx && return g.M[1, :]
    field === :Myy && return g.M[2, :]
    field === :Mxy && return g.M[3, :]
    field === :Qx && return g.Q[1, :]
    field === :Qy && return g.Q[2, :]
    error("unknown element field :$field")
end

"""
per-node values of `field` over the shell elements `ks` (columns of `g.quads`):
element fields are averaged ONLY over those elements (a wall's contour is not
polluted by the slab framing into it); `:w`/`:umag` are nodal displacements
"""
function nodal_field(g::ShellGeo, field::Symbol, ks = axes(g.quads, 2))
    nn = size(g.nodes, 2)
    if field === :umag
        return [norm(g.disp[:, i]) for i in 1:nn]
    elseif field === :w
        # out-of-plane displacement along each node's (averaged) element normal
        acc = zeros(nn); cnt = zeros(Int, nn)
        for k in ks, a in 1:4
            i = g.quads[a, k]
            acc[i] += dot(g.disp[:, i], g.normals[:, k]); cnt[i] += 1
        end
        return [c > 0 ? s / c : 0.0 for (s, c) in zip(acc, cnt)]
    end
    vals = element_field(g, field)
    acc = zeros(nn); cnt = zeros(Int, nn)
    for k in ks, a in 1:4
        acc[g.quads[a, k]] += vals[k]; cnt[g.quads[a, k]] += 1
    end
    return [c > 0 ? s / c : 0.0 for (s, c) in zip(acc, cnt)]
end

# columns of g.quads belonging to a set of elements
shell_columns(g::ShellGeo, model, els) = [findfirst(==(findfirst(e -> e === el, model.elements)), g.elements) for el in els]

# quads → triangles for mesh!
tri_faces(q::AbstractMatrix{Int}, ks) =
    reduce(vcat, [[q[1, k] q[2, k] q[3, k]; q[1, k] q[3, k] q[4, k]] for k in ks])

# quad edges as line segments
quad_edges(P, q, ks) = reduce(vcat, [P[q[[a, mod1(a + 1, 4)], k]] for k in ks for a in 1:4])

# ── wall elevation / slab plan ───────────────────────────────────────────────

# wall-plane coordinates (s, z) of global points (columns of X)
wall_sz(w::WallMesh, X::AbstractMatrix) =
    [Point2f(dot(w.along, X[:, i] - w.origin), dot(w.up, X[:, i] - w.origin)) for i in axes(X, 2)]

# a contoured shell drawing: mesh coloured by `vals`, element edges, colourbar beside it
function _contour_figure(P0, P, quads, ks, vals, field; title, colormap, colorrange, width, maxheight)
    # keep only this drawing's nodes: the geometry holds every node of the model, and mesh!
    # (hence the limits) would take them all in
    used = sort!(unique(vec(quads[:, ks])))
    remap = Dict(i => k for (k, i) in enumerate(used))
    quads = [remap[quads[a, k]] for a in 1:4, k in ks]
    ks = axes(quads, 2)
    P = P[used]; vals = vals[used]; P0 = P0 === nothing ? nothing : P0[used]
    fig = Figure()
    ax = Axis(fig[1, 1]; title = title)
    cm, cr = field_colormap(field, vals)
    colormap === nothing || (cm = colormap)
    colorrange === nothing || (cr = colorrange)
    P0 === nothing || linesegments!(ax, quad_edges(P0, quads, ks); color = kjl_lightgray, linewidth = 0.4)
    m = mesh!(ax, P, tri_faces(quads, ks); color = vals, colormap = cm, colorrange = cr, shading = NoShading)
    linesegments!(ax, quad_edges(P, quads, ks); color = (kjl_black, 0.35), linewidth = 0.3)
    Colorbar(fig[1, 2], m; label = get(FIELD_LABELS, field, string(field)))
    structurestyle!(ax; scale = true, unit = "m")
    fitstructure!(fig, ax; width = width - 62, maxheight = maxheight, padpt = 4)
    return fig
end

"""
    fig_wall(model, w; field = :N2, scale = 0.0, title, geo, colormap, colorrange, width = W_PLAN)

Elevation of the wall `w` (a `WallMesh` inside `model`): the mesh deformed by
`scale × displacement` (in-plane components), coloured by `field` averaged to
the nodes over this wall's elements (forces red − / white 0 / blue +; displacements white →
orange); the undeformed mesh in light gray.
"""
function fig_wall(model, w; field = :N2, scale = 0.0, title = "", geo = ShellGeo(model),
    colormap = nothing, colorrange = nothing, width = W_PLAN, maxheight = 1.1width)
    ks = shell_columns(geo, model, w.elements)
    vals = field_scale(field) .* nodal_field(geo, field, ks)
    P0 = wall_sz(w, geo.nodes)
    P = wall_sz(w, geo.nodes .+ scale .* geo.disp)
    return _contour_figure(scale == 0 ? nothing : P0, P, geo.quads, ks, vals, field;
        title, colormap, colorrange, width, maxheight)
end

"""
    fig_slab(model; field = :Mxx, elements = all shells, title, geo, colormap, width = W_HALF)

Plan view (x, y) of horizontal shell elements, coloured by `field`.
"""
function fig_slab(model; field = :Mxx, elements = nothing, title = "", geo = ShellGeo(model),
    colormap = nothing, colorrange = nothing, width = W_HALF, maxheight = 1.1width)
    ks = elements === nothing ? collect(axes(geo.quads, 2)) : shell_columns(geo, model, elements)
    vals = field_scale(field) .* nodal_field(geo, field, ks)
    P = [Point2f(geo.nodes[1, i], geo.nodes[2, i]) for i in axes(geo.nodes, 2)]
    return _contour_figure(nothing, P, geo.quads, ks, vals, field; title, colormap, colorrange, width, maxheight)
end

"""
    fig_principal(model, w; geo, title, maxlen, width = W_PLAN)

Principal membrane-force trajectories of wall `w`: at each element centre a
segment along each principal direction, length ∝ √|N|, blue = tension (ties),
red = compression (struts). The strut-and-tie reading of a wall.
"""
function fig_principal(model, w; geo = ShellGeo(model), title = "principal membrane forces", maxlen = nothing,
    width = W_PLAN, maxheight = 1.1width)
    ks = shell_columns(geo, model, w.elements)
    Nmax = maximum(k -> max(abs(geo.N1[k]), abs(geo.N2[k])), ks)
    h = maxlen === nothing ? 0.9 * minimum(diff(w.s)) : maxlen
    comp = Point2f[]; tens = Point2f[]
    for k in ks
        c = geo.centers[:, k]
        R = geo.frames[k]
        for (Nk, φ) in ((geo.N1[k], geo.θ[k]), (geo.N2[k], geo.θ[k] + π / 2))
            d = cos(φ) * R[1, :] + sin(φ) * R[2, :]               # global direction in the shell plane
            ℓ = h * sqrt(abs(Nk) / Nmax) / 2                     # √-scaled so weak regions stay visible
            seg = wall_sz(w, hcat(c - ℓ * d, c + ℓ * d))
            append!(Nk < 0 ? comp : tens, seg)
        end
    end
    fig = Figure()
    ax = Axis(fig[1, 1]; title = title)
    P0 = wall_sz(w, geo.nodes)
    linesegments!(ax, quad_edges(P0, geo.quads, ks); color = (kjl_lightgray, 0.6), linewidth = 0.3)
    linesegments!(ax, tens; color = kjl_tension, linewidth = 1.3)
    linesegments!(ax, comp; color = kjl_compression, linewidth = 1.3)
    Legend(fig[1, 2], [LineElement(color = kjl_tension, linewidth = 2), LineElement(color = kjl_compression, linewidth = 2)],
        ["tension (ties)", "compression (struts)"])
    structurestyle!(ax; scale = true, unit = "m")
    fitstructure!(fig, ax; width = width - 95, maxheight = maxheight, padpt = 4)
    return fig
end

# ── bars and series ──────────────────────────────────────────────────────────

"""
    fig_bars(labels, values; ylabel, title, highlight, alert, rotation, width = W_HALF)

Bar chart. All bars kjl_blue, unless `highlight` (index or indices: the bar the figure is
about, kjl_emphasis) or `alert` (a bad bar, kjl_alert) is given; the others are then gray.
"""
function fig_bars(labels, values; ylabel = "", title = "", highlight = nothing, alert = nothing,
    rotation = 0.0, width = W_HALF, height = 0.95width)
    n = length(values)
    colors = if highlight === nothing && alert === nothing
        fill(kjl_emphasis, n)
    else
        c = fill(kjl_lightgray, n)
        highlight === nothing || (c[collect(highlight)] .= kjl_emphasis)
        alert === nothing || (c[collect(alert)] .= kjl_alert)
        c
    end
    fig = Figure(size = (width, height))
    ax = Axis(fig[1, 1]; ylabel = ylabel, title = title, xticks = (1:n, labels),
        xticklabelrotation = rotation, xticksvisible = false)
    barplot!(ax, 1:n, values; color = colors)
    hlines!(ax, [0.0]; color = kjl_black, linewidth = 0.6, linestyle = :solid)
    return fig
end

"""
    fig_series(xs, series; xlabel, ylabel, ref, title, xscale, width = W_HALF)

`series` = vector of `(label, ys)`; the FIRST series is the emphasised one (kjl_blue), the
rest grayscale with distinct linestyles/markers (`kjl_series`). `ref` an optional horizontal
reference value (gray dashed).
"""
function fig_series(xs, series; xlabel = "", ylabel = "", ref = nothing, title = "", xscale = identity,
    reflabel = "reference", legendpos = :rb, width = W_HALF, height = 0.95width)
    fig = Figure(size = (width, height))
    ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, title = title, xscale = xscale)
    ref === nothing || hlines!(ax, [ref]; label = reflabel)
    for (i, (lab, ys)) in enumerate(series)
        scatterlines!(ax, xs, ys; label = lab, kjl_series(i; marker = true)...)
    end
    (length(series) > 1 || ref !== nothing) && axislegend(ax; position = legendpos)
    return fig
end

# ── 3D view ──────────────────────────────────────────────────────────────────

"""
    fig_3d(model; field = :umag, scale = 0.0, title, azimuth, elevation, linecolor = :frame, width = W_PLAN)

Shells (coloured by `field`, deformed by `scale × u`) and line elements
(deformed by their end displacements) in one 3D view on a ground plane, true proportions;
the undeformed line elements in light gray. `linecolor = :axial` colours members blue
(tension) / red (compression) with width ∝ |N|; otherwise dark gray.
"""
function fig_3d(model; field = :umag, scale = 0.0, title = "", geo = ShellGeo(model),
    azimuth = 1.1π, elevation = 0.35, colormap = nothing, width = W_PLAN, ratio = 0.72, linecolor = :frame,
    nodes = true, halo = true)
    size_shells(geo) == 0 && return _fig_3d_lines(model; scale, title, azimuth, elevation, width, ratio,
        linecolor, nodes, halo)
    # shells + line elements: hidden-surface vector drawing (kjlMakie `hiddenface!`): faces
    # painter-sorted, members and nodes behind a nearer face removed, crossings as in hiddenline!
    X = geo.nodes .+ scale .* geo.disp
    P = [Point3d(X[:, i]...) for i in axes(X, 2)]
    ks = collect(axes(geo.quads, 2))
    vals = field_scale(field) .* nodal_field(geo, field, ks)
    used = unique(vec(geo.quads))
    cm, cr = field_colormap(field, vals[used])
    colormap === nothing || (cm = colormap)
    tris = tri_faces(geo.quads, ks)
    es = Point3d[]; ee = Point3d[]
    for k in ks, a in 1:4
        push!(es, P[geo.quads[a, k]]); push!(ee, P[geo.quads[mod1(a + 1, 4), k]])
    end
    lines_ = [el for el in model.elements if el isa Union{FrameElement,TrussElement}]
    S = [P[el.nodeStart.index] for el in lines_]; E = [P[el.nodeEnd.index] for el in lines_]
    if linecolor === :axial
        N = [axial_force(model.results, el) for el in lines_]
        Nmax = max(maximum(abs, N; init = 0.0), eps())
        C = force_color.(N); W = [0.4 + 2.6 * abs(n) / Nmax for n in N]
    else
        C = fill(kjl_darkgray, length(lines_)); W = fill(1.2, length(lines_))
    end
    nds = nodes ? unique(vcat(S, E)) : nothing
    fig = Figure()
    ax = Axis(fig[1, 1]; title = title)
    m = hiddenface!(ax, P, tris, S, E; azimuth, elevation, facecolor = vals, colormap = cm, colorrange = cr,
        edges = (es, ee), color = C, linewidth = W, nodes = nds, nodesize = 2.5, nodestrokecolor = kjl_black,
        nodestrokewidth = 0.5, halo, ground = true, width = width - 60)
    Colorbar(fig[1, 2], m; label = get(FIELD_LABELS, field, string(field)), height = Relative(0.6))
    structurestyle!(ax)
    fitstructure!(fig, ax; width = width - 68, maxheight = ratio * width, padpt = 4)
    return fig
end

size_shells(g::ShellGeo) = size(g.quads, 2)

# line-only models: hidden-line vector drawing (kjlMakie `hiddenline!`), correct occlusion in PDF
function _fig_3d_lines(model; scale, title, azimuth, elevation, width, ratio, linecolor, nodes, halo)
    lines_ = [el for el in model.elements if el isa Union{FrameElement,TrussElement}]
    u(nd) = scale == 0 ? zeros(3) : displacement(model.results, nd)[1:3]
    S = Point3d[]; E = Point3d[]; C = Any[]; W = Float64[]
    N = linecolor === :axial ? [axial_force(model.results, el) for el in lines_] : zeros(length(lines_))
    Nmax = max(maximum(abs, N; init = 0.0), eps())
    for (el, n) in zip(lines_, N)
        push!(S, Point3d((el.nodeStart.position .+ scale .* u(el.nodeStart))...))
        push!(E, Point3d((el.nodeEnd.position .+ scale .* u(el.nodeEnd))...))
        push!(C, linecolor === :axial ? force_color(n) : kjl_darkgray)
        push!(W, linecolor === :axial ? 0.4 + 2.6 * abs(n) / Nmax : 1.2)
    end
    used = unique([nd for el in lines_ for nd in (el.nodeStart, el.nodeEnd)])
    nds = nodes ? [Point3d((nd.position .+ scale .* u(nd))...) for nd in used] : nothing
    fig = Figure()
    ax = Axis(fig[1, 1]; title = title)
    if scale != 0                                              # undeformed reference, light, behind
        P0, _ = project3d([Point3d(x.position...) for el in lines_ for x in (el.nodeStart, el.nodeEnd)], azimuth, elevation)
        linesegments!(ax, P0; color = kjl_lightgray, linewidth = 0.4)
    end
    hiddenline!(ax, S, E; azimuth, elevation, color = C, linewidth = W, nodes = nds, nodesize = 2.5,
        nodestrokecolor = kjl_black, nodestrokewidth = 0.5,
        halo, ground = true, width = width)
    structurestyle!(ax)
    fitstructure!(fig, ax; width = width - 8, maxheight = ratio * width, padpt = 4)
    return fig
end

# ── line elements (frames, trusses) ─────────────────────────────────────────

view2d(X, plane) = plane === :XZ ? Point2f(X[1], X[3]) : plane === :XY ? Point2f(X[1], X[2]) : Point2f(X[2], X[3])

# deformed polyline of a line element: exact curve for frames, straight for trusses
function member_curve(model, el, scale; resolution = 20)
    if el isa FrameElement
        ed = ElementDisplacements(el, model; resolution = resolution)
        return [ed.basepositions[:, k] .+ scale .* ed.uglobal[:, k] for k in axes(ed.basepositions, 2)]
    else
        u1 = displacement(model.results, el.nodeStart)[1:3]; u2 = displacement(model.results, el.nodeEnd)[1:3]
        return [el.nodeStart.position .+ scale .* u1, el.nodeEnd.position .+ scale .* u2]
    end
end

"""
    fig_lines(model; plane = :XZ, scale = 0.0, color = :axial, diagram = nothing, dscale = nothing,
              title, labels = false, widths = nothing, width = W_WIDE, scalebar = false)

2D view (`plane` = :XZ, :XY or :YZ) of the model's frame and truss members:
undeformed in light gray, deformed by `scale × u` along each member's EXACT
deflected curve, coloured by axial force (blue tension, red compression)
with line width ∝ |N| (or `widths`, one per member [pt], e.g. ∝ area), nodes as white dots. `diagram = :Mz` / `:My`
overlays that bending-moment diagram on each frame member's TENSION side (scaled by `dscale`):
black outline, 50 % black fill (one internal force → black). `labels = true` prints each member's N.
"""
function fig_lines(model; plane = :XZ, scale = 0.0, color = :axial, diagram = nothing, dscale = nothing,
    title = "", labels = false, widths = nothing, width = W_WIDE, maxheight = 0.7width, scalebar = false,
    nodes = true)
    lines_ = [el for el in model.elements if el isa Union{FrameElement,TrussElement}]
    fig = Figure()
    ax = Axis(fig[1, 1]; title = title)
    for el in lines_
        lines!(ax, [view2d(el.nodeStart.position, plane), view2d(el.nodeEnd.position, plane)]; color = kjl_lightgray, linewidth = 0.6)
    end
    N = [axial_force(model.results, el) for el in lines_]
    Nmax = maximum(abs, N; init = 0.0); Nmax = Nmax == 0 ? 1.0 : Nmax
    if diagram !== nothing
        frames = [el for el in lines_ if el isa FrameElement]
        fs = [InternalForces(model, el; resolution = 30) for el in frames]
        Ms = [getproperty(f, diagram) for f in fs]
        Mmax = maximum(m -> maximum(abs, m), Ms; init = 0.0)
        span = maximum(el -> norm(el.nodeEnd.position - el.nodeStart.position), frames)
        ds = dscale === nothing ? 0.15 * span / max(Mmax, eps()) : dscale
        fills = Vector{Point2f}[]
        for (el, f, M) in zip(frames, fs, Ms)
            a, b = view2d(el.nodeStart.position, plane), view2d(el.nodeEnd.position, plane)
            # drawn on the TENSION side: positive (sagging) Mz is tension on the member's −y face
            # (−z for My), so offset along −(local y / z) projected into the view plane
            Λ = local_frame(el.nodeStart.position, el.nodeEnd.position, el.rollangle)
            ax_ = diagram === :Mz ? Λ[2, :] : Λ[3, :]
            nrm = -view2d(ax_, plane)
            norm(nrm) > 1e-6 || (t = (b - a) / norm(b - a); nrm = Point2f(t[2], -t[1]))
            nrm = nrm / norm(nrm)
            L = norm(el.nodeEnd.position - el.nodeStart.position)
            base = [a + (x / L) * (b - a) for x in f.x]
            curve = [p + ds * m * nrm for (p, m) in zip(base, M)]
            # one simple polygon per run of same-sign moment (split at zero crossings):
            # per-segment quads would show hairline seams in vector output
            run_b = [base[1]]; run_c = [curve[1]]
            for k in 2:length(base)
                if M[k-1] * M[k] < 0
                    t = M[k-1] / (M[k-1] - M[k]); z = base[k-1] + t * (base[k] - base[k-1])
                    push!(run_b, z); push!(run_c, z)
                    push!(fills, vcat(run_b, reverse(run_c)))
                    run_b = [z]; run_c = [z]
                end
                push!(run_b, base[k]); push!(run_c, curve[k])
            end
            push!(fills, vcat(run_b, reverse(run_c)))
            lines!(ax, curve; color = kjl_black, linewidth = 0.8)
        end
        isempty(fills) || poly!(ax, fills; color = (kjl_black, 0.5), strokewidth = 0)
    end
    for (k, (el, n)) in enumerate(zip(lines_, N))
        pts = [view2d(x, plane) for x in member_curve(model, el, scale)]
        c = color === :axial ? force_color(n) : color
        lw = widths !== nothing ? widths[k] : color === :axial ? 0.75 + 3.5 * abs(n) / Nmax : 1.5
        lines!(ax, pts; color = c, linewidth = lw)
    end
    if nodes
        nds = unique([nd for el in lines_ for nd in (el.nodeStart, el.nodeEnd)])
        u(nd) = displacement(model.results, nd)[1:3]
        scatter!(ax, [view2d(nd.position .+ scale .* u(nd), plane) for nd in nds];
            color = :white, strokecolor = kjl_black, strokewidth = 0.6, markersize = 4)
    end
    if labels
        for el in lines_
            n = round(axial_force(model.results, el), digits = 1)
            text!(ax, view2d((el.nodeStart.position + el.nodeEnd.position) / 2, plane);
                text = replace(string(n), "-" => "−"), fontsize = 6.5, color = kjl_darkgray, align = (:center, :bottom))
        end
    end
    color === :axial && Legend(fig[2, 1], [LineElement(color = kjl_tension, linewidth = 2.5), LineElement(color = kjl_compression, linewidth = 2.5)],
        ["tension", "compression"]; orientation = :horizontal, tellwidth = false)
    structurestyle!(ax; scale = scalebar, unit = "m")
    fitstructure!(fig, ax; width = width - 8, maxheight = maxheight, padpt = labels ? 14 : 6)
    return fig
end

"""
    fig_member(model, el; title, width = W_HALF)

Axial force, shear and bending-moment diagrams along one member (exact
recovery: `InternalForces`, stations at load breakpoints), stacked; each a black outline with
a 50 % black fill (`diagramband!`).
"""
function fig_member(model, el; title = "", width = W_HALF, height = 1.2width)
    f = InternalForces(model, el; resolution = 60)
    fig = Figure(size = (width, height))
    axs = Axis[]
    for (row, (y, lab)) in enumerate(((f.N, "N  [kN]"), (f.Vy, rich("V", subscript("y"), "  [kN]")),
                                      (f.Mz, rich("M", subscript("z"), "  [kN·m]"))))
        ax = Axis(fig[row, 1]; aspect = nothing, ylabel = lab, xlabel = row == 3 ? "{x} along member [m]" : "",
            title = row == 1 ? title : "", xticklabelsvisible = row == 3)
        diagramband!(ax, f.x, y)
        push!(axs, ax)
    end
    linkxaxes!(axs...)
    rowgap!(fig.layout, 6)
    return fig
end

println("figures.jl loaded")
