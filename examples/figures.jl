#=
Figure generators for the shell walkthroughs.

Plain functions of Asap's plot-extraction objects (`ShellGeo`, `WallMesh`,
the solved `Model`) that return Makie figures in the kjlMakie house style:

    fig_wall(model, w; field = :N2, scale = 0.0)    wall elevation, contoured field on the (scaled) deformed mesh
    fig_principal(model, w)                         principal membrane-force trajectories (struts / ties)
    fig_bars(labels, values; ylabel)                bar chart (pier shears, storey shears, …)
    fig_slab(model; field = :Mxx, elements)         plan view of a slab, contoured
    fig_3d(model; field = :N2, scale)               3D view of shells + frames, deformed
    fig_series(xs, series; xlabel, ylabel, ref)     convergence / parametric curves
    fig_lines(model; plane, scale, diagram)         2D view of frames/trusses: deformed (exact member curves),
                                                    coloured by axial force, optional moment diagram
    fig_member(model, el)                           N / V / M diagrams along one member

`include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)` from any
walkthrough (it loads the plotting stack); figures go to `figures/` next to the walkthrough.
`report_data(name; key = value, …)` writes the numbers a walkthrough's Typst report prints
to `reports/data/<name>.json`.
=#

using Asap, LinearAlgebra
using kjlMakie, CairoMakie
CairoMakie.activate!()
set_theme!(kjl_light)
update_theme!(backgroundcolor = :white)        # PNGs: opaque white, not transparent

# each walkthrough saves into `figures/` next to itself: call `set_figdir(@__DIR__)`
const FIGDIR = Ref(joinpath(@__DIR__, "figures"))
set_figdir(dir) = (FIGDIR[] = joinpath(dir, "figures"); mkpath(FIGDIR[]); FIGDIR[])
# Saved PNGs are cropped to their content plus a margin: axes with a data aspect (elevations,
# plans) otherwise leave wide empty bands, which waste page space in the reports.
function savefig(name, fig; pad = 24)
    path = joinpath(FIGDIR[], name)
    save(path, fig; px_per_unit = 2)
    endswith(name, ".png") && _crop_whitespace!(path; pad)
    return fig
end
function _crop_whitespace!(path; pad = 24)
    img = Makie.FileIO.load(path)
    C = Makie.Colors
    ink = map(c -> min(C.red(c), C.green(c), C.blue(c)) < 0.98, img)
    rows, cols = findall(vec(any(ink; dims = 2))), findall(vec(any(ink; dims = 1)))
    (isempty(rows) || isempty(cols)) && return path
    r = max(first(rows) - pad, 1):min(last(rows) + pad, size(img, 1))
    c = max(first(cols) - pad, 1):min(last(cols) + pad, size(img, 2))
    (length(r), length(c)) == size(img) || Makie.FileIO.save(path, img[r, c])
    return path
end

# ── report data ──────────────────────────────────────────────────────────────
# Each walkthrough ends with `report_data("A1"; key = value, …)`: the numbers its
# Typst report (`reports/`) prints, written to `reports/data/<name>.json`, so no
# number in a report is typed by hand. Values: numbers, strings, symbols,
# vectors and (named) tuples of those; non-finite floats are written as null.
using JSON
const REPORT_DATA = joinpath(@__DIR__, "reports", "data")
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

const FIELD_LABELS = Dict(
    :Nxx => "Nxx  [force/length]", :Nyy => "Nyy  [force/length]", :Nxy => "Nxy  [force/length]",
    :N1 => "N1 (major principal)  [force/length]", :N2 => "N2 (minor principal)  [force/length]",
    :Mxx => "Mxx  [force·length/length]", :Myy => "Myy  [force·length/length]",
    :Mxy => "Mxy  [force·length/length]", :Qx => "Qx  [force/length]", :Qy => "Qy  [force/length]",
    :w => "out-of-plane displacement  [length]", :umag => "|u|  [length]")

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

symmetric_range(v) = (m = maximum(abs, v; init = 0.0); m == 0 ? (-1.0, 1.0) : (-m, m))

# ── wall elevation ───────────────────────────────────────────────────────────

# figure size following the wall's aspect ratio (plot area ≈ 420 px tall + colourbar/legend)
wall_figsize(w::WallMesh) = (round(Int, clamp(420 * w.s[end] / w.z[end], 160, 900)) + 260, 480)

# limits from the plotted points only (the geometry holds every node of the
# model, so auto-limits would take in other walls), and the axis column sized
# to the data aspect so a colourbar/legend sits right beside the wall
function fit_wall_axis!(fig, ax, pts; pad = 0.03)
    s = extrema(p -> p[1], pts); z = extrema(p -> p[2], pts)
    ds, dz = s[2] - s[1], z[2] - z[1]
    m = pad * max(ds, dz)
    limits!(ax, s[1] - m, s[2] + m, z[1] - m, z[2] + m)
    colsize!(fig.layout, 1, Aspect(1, (ds + 2m) / (dz + 2m)))
    resize_to_layout!(fig)
    return fig
end

# wall-plane coordinates (s, z) of global points (columns of X)
wall_sz(w::WallMesh, X::AbstractMatrix) =
    [Point2f(dot(w.along, X[:, i] - w.origin), dot(w.up, X[:, i] - w.origin)) for i in axes(X, 2)]

"""
    fig_wall(model, w; field = :N2, scale = 0.0, title, geo, colormap)

Elevation of the wall `w` (a `WallMesh` inside `model`): the mesh deformed by
`scale × displacement` (in-plane components), coloured by `field` averaged to
the nodes over this wall's elements; the undeformed outline in grey.
"""
function fig_wall(model, w; field = :N2, scale = 0.0, title = "", geo = ShellGeo(model),
    colormap = :balance, colorrange = nothing)
    ks = shell_columns(geo, model, w.elements)
    vals = nodal_field(geo, field, ks)
    P0 = wall_sz(w, geo.nodes)
    P = wall_sz(w, geo.nodes .+ scale .* geo.disp)
    F = tri_faces(geo.quads, ks)
    fig = Figure(size = wall_figsize(w))
    ax = Axis(fig[1, 1]; aspect = DataAspect(), title = title, xlabel = "s  [length]", ylabel = "z  [length]")
    edges0 = [P0[geo.quads[[a, mod1(a + 1, 4)], k]] for k in ks for a in 1:4]
    linesegments!(ax, reduce(vcat, edges0); color = (:gray, 0.4), linewidth = 0.5)
    cr = colorrange === nothing ? (field in (:umag,) ? (0.0, maximum(vals[unique(vec(geo.quads[:, ks]))])) :
                                   symmetric_range(vals[unique(vec(geo.quads[:, ks]))])) : colorrange
    m = mesh!(ax, P, F; color = vals, colormap = colormap, colorrange = cr, shading = NoShading)
    edges = [P[geo.quads[[a, mod1(a + 1, 4)], k]] for k in ks for a in 1:4]
    linesegments!(ax, reduce(vcat, edges); color = :black, linewidth = 0.4)
    Colorbar(fig[1, 2], m; label = get(FIELD_LABELS, field, string(field)))
    used = unique(vec(geo.quads[:, ks]))
    return fit_wall_axis!(fig, ax, vcat(P0[used], P[used]))
end

"""
    fig_principal(model, w; geo, title, maxlen)

Principal membrane-force trajectories of wall `w`: at each element centre a
segment along each principal direction, length ∝ |N|, blue = compression
(struts), pink = tension (ties). The strut-and-tie reading of a wall.
"""
function fig_principal(model, w; geo = ShellGeo(model), title = "principal membrane forces", maxlen = nothing)
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
    fig = Figure(size = wall_figsize(w))
    ax = Axis(fig[1, 1]; aspect = DataAspect(), title = title, xlabel = "s  [length]", ylabel = "z  [length]")
    P0 = wall_sz(w, geo.nodes)
    linesegments!(ax, reduce(vcat, [P0[geo.quads[[a, mod1(a + 1, 4)], k]] for k in ks for a in 1:4]);
        color = (:gray, 0.35), linewidth = 0.5)
    linesegments!(ax, comp; color = kjl_blue, linewidth = 1.6, label = "compression")
    linesegments!(ax, tens; color = kjl_pink, linewidth = 1.6, label = "tension")
    Legend(fig[1, 2], ax; framevisible = false)
    return fit_wall_axis!(fig, ax, P0[unique(vec(geo.quads[:, ks]))])
end

# ── bars and series ──────────────────────────────────────────────────────────

function fig_bars(labels, values; ylabel = "", title = "", color = kjl_blue, rotation = 0.0)
    fig = Figure(size = (max(520, 46 * length(labels)), 360))
    ax = Axis(fig[1, 1]; ylabel = ylabel, title = title, xticks = (1:length(labels), labels),
        xticklabelrotation = rotation)
    barplot!(ax, 1:length(values), values; color = color, gap = 0.25)
    hlines!(ax, [0.0]; color = :black, linewidth = 0.5)
    return fig
end

"""
    fig_series(xs, series; xlabel, ylabel, ref, title, xscale)

`series` = vector of `(label, ys)`; `ref` an optional horizontal reference value.
"""
function fig_series(xs, series; xlabel = "", ylabel = "", ref = nothing, title = "", xscale = identity,
    reflabel = "reference", legendpos = :rb)
    fig = Figure(size = (560, 380))
    ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, title = title, xscale = xscale)
    ref === nothing || hlines!(ax, [ref]; color = :black, linestyle = :dash, label = reflabel)
    colors = (kjl_blue, kjl_pink, kjl_green, kjl_orange, kjl_turquoise)
    for (i, (lab, ys)) in enumerate(series)
        scatterlines!(ax, xs, ys; label = lab, color = colors[mod1(i, length(colors))])
    end
    axislegend(ax; position = legendpos)
    return fig
end

# ── slab plan view ───────────────────────────────────────────────────────────

"""
    fig_slab(model; field = :Mxx, elements = all shells, title)

Plan view (x, y) of horizontal shell elements, coloured by `field`.
"""
function fig_slab(model; field = :Mxx, elements = nothing, title = "", geo = ShellGeo(model), colormap = :balance)
    ks = elements === nothing ? collect(axes(geo.quads, 2)) : shell_columns(geo, model, elements)
    vals = nodal_field(geo, field, ks)
    P = [Point2f(geo.nodes[1, i], geo.nodes[2, i]) for i in axes(geo.nodes, 2)]
    F = tri_faces(geo.quads, ks)
    fig = Figure(size = (560, 480))
    ax = Axis(fig[1, 1]; aspect = DataAspect(), title = title, xlabel = "x  [length]", ylabel = "y  [length]")
    used = unique(vec(geo.quads[:, ks]))
    cr = field === :umag ? (0.0, maximum(vals[used])) : symmetric_range(vals[used])
    m = mesh!(ax, P, F; color = vals, colormap = colormap, colorrange = cr, shading = NoShading)
    linesegments!(ax, reduce(vcat, [P[geo.quads[[a, mod1(a + 1, 4)], k]] for k in ks for a in 1:4]);
        color = (:black, 0.3), linewidth = 0.4)
    Colorbar(fig[1, 2], m; label = get(FIELD_LABELS, field, string(field)))
    return fit_wall_axis!(fig, ax, P[used])
end

# ── 3D view ──────────────────────────────────────────────────────────────────

"""
    fig_3d(model; field = :umag, scale = 0.0, title, azimuth, elevation)

Shells (coloured by `field`, deformed by `scale × u`) and line elements
(dark lines, deformed by their end displacements) in one 3D view; the
undeformed line elements in grey.
"""
function fig_3d(model; field = :umag, scale = 0.0, title = "", geo = ShellGeo(model),
    azimuth = 1.1π, elevation = 0.35, colormap = :viridis, size = (760, 640), linecolor = :frame)
    fig = Figure(size = size)
    ax = Axis3(fig[1, 1]; aspect = :data, title = title, azimuth = azimuth, elevation = elevation)
    X0 = geo.nodes
    X = geo.nodes .+ scale .* geo.disp
    P = [Point3f(X[:, i]) for i in axes(X, 2)]
    if size_shells(geo) > 0
        ks = axes(geo.quads, 2)
        vals = nodal_field(geo, field, ks)
        used = unique(vec(geo.quads))
        cr = field === :umag ? (0.0, max(maximum(vals[used]), eps())) : symmetric_range(vals[used])
        m = mesh!(ax, P, tri_faces(geo.quads, ks); color = vals, colormap = colormap, colorrange = cr)
        linesegments!(ax, reduce(vcat, [P[geo.quads[[a, mod1(a + 1, 4)], k]] for k in ks for a in 1:4]);
            color = (:black, 0.25), linewidth = 0.4)
        Colorbar(fig[1, 2], m; label = get(FIELD_LABELS, field, string(field)))
    end
    lines_ = [el for el in model.elements if el isa Union{FrameElement,TrussElement}]
    if !isempty(lines_)
        seg0 = Point3f[]; seg = Point3f[]
        for el in lines_
            i, j = el.nodeStart.index, el.nodeEnd.index
            push!(seg0, Point3f(X0[:, i]), Point3f(X0[:, j]))
            push!(seg, P[i], P[j])
        end
        scale == 0 || linesegments!(ax, seg0; color = (:gray, 0.5), linewidth = 1)
        if linecolor === :axial
            N = [axial_force(model.results, el) for el in lines_]
            Nmax = max(maximum(abs, N), eps())
            cols = reduce(vcat, [[n < 0 ? kjl_blue : kjl_pink, n < 0 ? kjl_blue : kjl_pink] for n in N])
            widths = reduce(vcat, [[0.5 + 4 * abs(n) / Nmax, 0.5 + 4 * abs(n) / Nmax] for n in N])
            linesegments!(ax, seg; color = cols, linewidth = widths)
        else
            linesegments!(ax, seg; color = kjl_darkgray, linewidth = 2.5)
        end
    end
    return fig
end

size_shells(g::ShellGeo) = size(g.quads, 2)

# ── line elements (frames, trusses) ─────────────────────────────────────────

view2d(X, plane) = plane === :XZ ? Point2f(X[1], X[3]) : plane === :XY ? Point2f(X[1], X[2]) : Point2f(X[2], X[3])
const PLANE_LABELS = Dict(:XZ => ("x  [length]", "z  [length]"), :XY => ("x  [length]", "y  [length]"),
    :YZ => ("y  [length]", "z  [length]"))

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
    fig_lines(model; plane = :XZ, scale = 0.0, color = :axial, diagram = nothing, dscale = nothing, title)

2D view (`plane` = :XZ, :XY or :YZ) of the model's frame and truss members:
undeformed in grey, deformed by `scale × u` along each member's EXACT
deflected curve, coloured by axial force (blue compression, pink tension)
with line width ∝ |N| (or `widths`, one per member, e.g. ∝ area). `diagram = :Mz` / `:My` overlays that bending-moment
diagram normal to each frame member (scaled by `dscale`).
"""
function fig_lines(model; plane = :XZ, scale = 0.0, color = :axial, diagram = nothing, dscale = nothing,
    title = "", size = (760, 460), labels = false, widths = nothing)
    lines_ = [el for el in model.elements if el isa Union{FrameElement,TrussElement}]
    fig = Figure(size = size)
    xl, yl = PLANE_LABELS[plane]
    ax = Axis(fig[1, 1]; aspect = DataAspect(), title = title, xlabel = xl, ylabel = yl,
        xautolimitmargin = (0.06, 0.06), yautolimitmargin = (0.15, 0.15))
    for el in lines_
        lines!(ax, [view2d(el.nodeStart.position, plane), view2d(el.nodeEnd.position, plane)]; color = (:gray, 0.45), linewidth = 1)
    end
    N = [axial_force(model.results, el) for el in lines_]
    Nmax = maximum(abs, N; init = 0.0); Nmax = Nmax == 0 ? 1.0 : Nmax
    for (k, (el, n)) in enumerate(zip(lines_, N))
        pts = [view2d(x, plane) for x in member_curve(model, el, scale)]
        c = color === :axial ? (n < 0 ? kjl_blue : kjl_pink) : color
        lw = widths !== nothing ? widths[k] : color === :axial ? 1 + 5 * abs(n) / Nmax : 2
        lines!(ax, pts; color = c, linewidth = lw)
    end
    if diagram !== nothing
        frames = [el for el in lines_ if el isa FrameElement]
        Ms = [getproperty(InternalForces(model, el; resolution = 30), diagram) for el in frames]
        Mmax = maximum(m -> maximum(abs, m), Ms; init = 0.0)
        span = maximum(el -> norm(el.nodeEnd.position - el.nodeStart.position), frames)
        ds = dscale === nothing ? 0.25 * span / max(Mmax, eps()) : dscale
        for (el, M) in zip(frames, Ms)
            f = InternalForces(model, el; resolution = 30)
            a, b = view2d(el.nodeStart.position, plane), view2d(el.nodeEnd.position, plane)
            t = (b - a) / norm(b - a); nrm = Point2f(-t[2], t[1])
            L = norm(el.nodeEnd.position - el.nodeStart.position)
            base = [a + (x / L) * (b - a) for x in f.x]
            curve = [p + ds * m * nrm for (p, m) in zip(base, M)]
            poly!(ax, vcat(base, reverse(curve)); color = (kjl_green, 0.25))
            lines!(ax, curve; color = kjl_green, linewidth = 1.2)
        end
    end
    if labels
        for el in lines_
            text!(ax, view2d((el.nodeStart.position + el.nodeEnd.position) / 2, plane);
                text = string(round(axial_force(model.results, el), digits = 1)), fontsize = 9, align = (:center, :bottom))
        end
    end
    color === :axial && Legend(fig[1, 2], [LineElement(color = kjl_blue, linewidth = 3), LineElement(color = kjl_pink, linewidth = 3)] .|> identity,
        ["compression", "tension"]; framevisible = false)
    resize_to_layout!(fig)
    return fig
end

"""
    fig_member(model, el; title)

Axial force, shear and bending-moment diagrams along one member (exact
recovery: `InternalForces`, stations at load breakpoints).
"""
function fig_member(model, el; title = "")
    f = InternalForces(model, el; resolution = 60)
    fig = Figure(size = (620, 560))
    for (row, (y, lab, c)) in enumerate(((f.N, "N  [force]", kjl_darkgray), (f.Vy, "Vy  [force]", kjl_turquoise),
                                          (f.Mz, "Mz  [force·length]", kjl_green)))
        ax = Axis(fig[row, 1]; ylabel = lab, xlabel = row == 3 ? "x along member  [length]" : "",
            title = row == 1 ? title : "")
        band!(ax, f.x, zeros(length(y)), y; color = (c, 0.25))
        lines!(ax, f.x, y; color = c, linewidth = 2)
        hlines!(ax, [0.0]; color = :black, linewidth = 0.5)
    end
    return fig
end

println("figures.jl loaded")
