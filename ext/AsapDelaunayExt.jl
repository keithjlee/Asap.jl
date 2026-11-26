"""
    AsapDelaunayExt

`mesh_polygon`: planar polygons with polygonal holes → triangular shell
meshes, through DelaunayTriangulation.jl (pure Julia: constrained Delaunay
triangulation with Ruppert-type quality refinement). Loaded when
DelaunayTriangulation is in the environment; core Asap stays dependency-free.

The polygons enter as a planar straight-line graph of CONSTRAINED SEGMENTS
(coincident vertices merged, edges split where another polygon's vertex lies
on them), not as boundary curves. That is what lets holes touch the outline
or each other — at a vertex, a vertex on an edge, or along a shared edge (a
notch drawn as a hole, adjacent openings drawn as two) — which boundary curves
cannot express. Interior `segments` (cut lines, load lines) enter the same
graph. After refinement the triangles outside the material (outside the
outline or inside a hole) are dropped, and the nodes on every input edge are
recovered geometrically.
"""
module AsapDelaunayExt

using Asap
using Asap: ShellMesh, AbstractShellSection
using LinearAlgebra: cross, dot
using DelaunayTriangulation
const DT = DelaunayTriangulation

# 2D points of a polygon as Float64 tuples, no closing repeat
function _poly(P)
    pts = P isa AbstractMatrix ? [(Float64(P[1, k]), Float64(P[2, k])) for k in axes(P, 2)] :
          [(Float64(p[1]), Float64(p[2])) for p in P]
    length(pts) > 1 && pts[end] == pts[1] && pop!(pts)
    length(pts) >= 3 || throw(ArgumentError("a polygon needs at least 3 vertices"))
    return pts
end

_area2(p) = sum(k -> p[k][1] * p[mod1(k + 1, length(p))][2] - p[mod1(k + 1, length(p))][1] * p[k][2], eachindex(p))

# even–odd point-in-polygon
function _inside(P, x, y)
    c = false
    n = length(P)
    for k in 1:n
        xi, yi = P[k]
        xj, yj = P[mod1(k + 1, n)]
        ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) && (c = !c)
    end
    return c
end

_cross(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])

# do segments ab and cd cross at a point interior to both? (touching at an
# endpoint and collinear overlap are not crossings)
function _crosses(a, b, c, d, tol)
    d1, d2 = _cross(a, b, c), _cross(a, b, d)
    d3, d4 = _cross(c, d, a), _cross(c, d, b)
    return ((d1 > tol && d2 < -tol) || (d1 < -tol && d2 > tol)) &&
           ((d3 > tol && d4 < -tol) || (d3 < -tol && d4 > tol))
end

# is p on edge k of poly (within tol)?
function _on_edge(p, poly, k, tol)
    a, b = poly[k], poly[mod1(k + 1, length(poly))]
    L = hypot(b[1] - a[1], b[2] - a[2])
    s = ((p[1] - a[1]) * (b[1] - a[1]) + (p[2] - a[2]) * (b[2] - a[2])) / L
    return abs(_cross(a, b, p)) / L <= tol && -tol <= s <= L + tol
end

function Asap.mesh_polygon(outline, section::AbstractShellSection{T};
    holes = [], segments = [], h = nothing, max_area = nothing, min_angle::Real = 25.0, points = [],
    origin = zeros(3), e1 = [1.0, 0, 0], e2 = [0.0, 1, 0], fixity::Symbol = :free,
    id::Symbol = :panel, seed::Integer = 1, reuse::AbstractVector = Asap.Node{T}[]) where {T}
    outer = _poly(outline)
    _area2(outer) < 0 && reverse!(outer)                    # outline counter-clockwise
    hs = map(holes) do H
        p = _poly(H)
        _area2(p) > 0 ? reverse(p) : p                       # holes clockwise
    end
    polys = (outer, hs...)
    scale = maximum(p -> maximum(q -> max(abs(q[1]), abs(q[2])), p), polys)
    tol = 1e-9 * max(scale, 1.0)
    for (k, H) in enumerate(hs)
        all(p -> _inside(outer, p...) || any(e -> _on_edge(p, outer, e, tol), eachindex(outer)), H) ||
            throw(ArgumentError("mesh_polygon: hole $k is not inside the outline"))
    end
    # vertices, coincident ones (touching polygons) merged
    pts = Tuple{Float64,Float64}[]
    function vid(p)
        i = findfirst(q -> abs(q[1] - p[1]) <= tol && abs(q[2] - p[2]) <= tol, pts)
        i === nothing || return i
        push!(pts, p)
        return length(pts)
    end
    ids = [[vid(p) for p in poly] for poly in polys]
    # in the material: inside or on the outline, not strictly inside a hole
    inmaterial(q) = (_inside(outer, q...) || any(k -> _on_edge(q, outer, k, tol), eachindex(outer))) &&
                    !any(H -> _inside(H, q...) && !any(k -> _on_edge(q, H, k, tol), eachindex(H)), hs)
    # extra points (point loads, ends of loaded lengths): must lie in the material
    for (k, p) in enumerate(points)
        q = (Float64(p[1]), Float64(p[2]))
        inmaterial(q) || throw(ArgumentError("mesh_polygon: point $k $(q) does not lie in the material"))
        vid(q)
    end
    # interior segments (cut lines, load lines): both ends in the material
    segs = map(enumerate(segments)) do (k, sg)
        a, b = (Float64(sg[1][1]), Float64(sg[1][2])), (Float64(sg[2][1]), Float64(sg[2][2]))
        a == b && throw(ArgumentError("mesh_polygon: segment $k has zero length"))
        inmaterial(a) && inmaterial(b) && inmaterial(((a[1] + b[1]) / 2, (a[2] + b[2]) / 2)) ||
            throw(ArgumentError("mesh_polygon: segment $k does not lie in the material"))
        (vid(a), vid(b))
    end
    # reused nodes in the panel's plane and in the material become vertices
    o, a, b = Asap._plane_axes(T, origin, e1, e2)
    nrm = cross(a, b)
    for n in reuse
        d = n.position - o
        abs(dot(d, nrm)) <= tol || continue
        q = (Float64(dot(d, a)), Float64(dot(d, b)))
        inmaterial(q) && vid(q)
    end
    # every input edge, split where another vertex lies on it (T-junctions); an
    # edge shared by two polygons (a hole along the outline or another hole) once
    edges = Set{NTuple{2,Int}}()
    raw = [[(poly[k], poly[mod1(k + 1, length(poly))]) for poly in ids for k in eachindex(poly)]; segs]
    chains = Vector{Int}[]             # each input edge with the vertices on it, in order
    for (i, j) in raw
        a, b = pts[i], pts[j]
        dx, dy = b[1] - a[1], b[2] - a[2]
        L = hypot(dx, dy)
        on = Tuple{Float64,Int}[]
        for (m, p) in enumerate(pts)
            (m == i || m == j) && continue
            s = ((p[1] - a[1]) * dx + (p[2] - a[2]) * dy) / L
            abs((p[1] - a[1]) * dy - (p[2] - a[2]) * dx) / L <= tol && tol < s < L - tol && push!(on, (s, m))
        end
        chain = [i; last.(sort!(on)); j]
        push!(chains, chain)
        for q in 1:length(chain)-1
            push!(edges, minmax(chain[q], chain[q+1]))
        end
    end
    # crossing segments: a hole crosses the outline or two holes overlap
    E = collect(edges)
    for x in eachindex(E), y in x+1:length(E)
        (a, b), (c, d) = E[x], E[y]
        _crosses(pts[a], pts[b], pts[c], pts[d], tol * scale) && throw(ArgumentError(
            "mesh_polygon: edges cross near $(pts[a]) → $(pts[b]): holes must lie inside the outline " *
            "and must not overlap each other (touching is allowed); segments must not cross an edge or " *
            "each other (split them at the crossing point)"))
    end
    tri = _triangulate(pts, ids[1], chains, tol, seed)
    # size: a number (uniform target edge length) or a function (x, y) -> h of the
    # in-plane position (graded meshes), both as the area of an equilateral triangle
    # of edge h; `max_area` caps either
    hfun = h isa Function
    amax = max_area !== nothing ? Float64(max_area) :
           h !== nothing && !hfun ? sqrt(3) / 4 * Float64(h)^2 : nothing
    # a cap on the vertex count turns any remaining runaway into an error; graded:
    # the expected count ∫ dA / (√3/4·h²) from a sampling of the material
    Amat = (_area2(outer) + sum(_area2, hs; init = 0.0)) / 2
    est = amax === nothing ? 0.0 : Amat / amax
    hfun && (est = max(est, _graded_count(h, outer, hs)))
    cap = max(20_000, 10 * length(pts) + ceil(Int, 20 * est))
    kw = (; min_angle = Float64(min_angle), max_points = cap, rng = DT.Random.Xoshiro(seed))
    amax === nothing || (kw = (; kw..., max_area = amax))
    if hfun
        # refine a material triangle whose area exceeds √3/4·h(centroid)²
        # (triangles inside holes are left alone: they are dropped)
        too_big = function (t, T)
            i, j, k = DT.triangle_vertices(T)
            p, q, r = DT.get_point(t, i, j, k)
            cx, cy = (p[1] + q[1] + r[1]) / 3, (p[2] + q[2] + r[2]) / 3
            any(H -> _inside(H, cx, cy), hs) && return false
            hc = Float64(h(cx, cy))
            hc > 0 || throw(ArgumentError("mesh_polygon: h($cx, $cy) = $hc must be positive"))
            A = abs((q[1] - p[1]) * (r[2] - p[2]) - (r[1] - p[1]) * (q[2] - p[2])) / 2
            return A > sqrt(3) / 4 * hc^2
        end
        kw = (; kw..., custom_constraint = too_big)
    end
    refine!(tri; kw...)
    DT.num_solid_vertices(tri) >= cap && throw(ErrorException(
        "mesh_polygon: refinement stopped at $cap vertices without meeting the quality targets. " *
        "Usual causes: two input vertices closer together than the mesh size (a point a hair off an " *
        "edge, near-duplicate polygon vertices), or a tiny `min_angle`-violating input feature."))
    # keep the material: inside the outline, outside every hole (centroid test)
    P = DT.get_points(tri)
    cells = Vector{Int}[]
    for t in DT.each_solid_triangle(tri)
        i, j, k = DT.triangle_vertices(t)
        cx, cy = (P[i][1] + P[j][1] + P[k][1]) / 3, (P[i][2] + P[j][2] + P[k][2]) / 3
        _inside(outer, cx, cy) && !any(H -> _inside(H, cx, cy), hs) && push!(cells, [i, j, k])
    end
    mesh = Asap.shell_mesh([P[i] for i in eachindex(P)], cells, section;
        origin = origin, e1 = e1, e2 = e2, fixity = fixity, id = id, reuse = reuse)
    # the mesh nodes on every input edge, in order along it
    boundary = [[_edge_nodes(mesh, poly[k], poly[mod1(k + 1, length(poly))], tol) for k in eachindex(poly)]
                for poly in polys]
    return ShellMesh(mesh.nodes, mesh.elements, mesh.points, boundary, mesh.origin, mesh.e1, mesh.e2)
end

# Triangulate the input graph (seeded: the same input gives the same mesh).
#
# The outline is the domain's boundary curve, so nothing outside it is
# triangulated or refined: with segments alone the domain is the convex hull,
# and a vertex a rounding error outside an edge (a load point on a sloped edge
# lands ±1e-16 off it) leaves a sliver between hull and outline that Ruppert
# refinement splits without end. Holes and interior segments stay constrained
# segments, which is what lets holes touch the outline and each other.
#
# Vertices lying ON an input edge (T-junctions: load points, touching hole
# corners, segment ends, reused nodes) are not triangulated as free points:
# off-line by rounding, a run of them would give zero-area slivers. They are
# skipped at first and then inserted by splitting the edge they lie on, as the
# refinement itself splits segments, which keeps them on the edge
# topologically whatever their rounding.
function _triangulate(pts, outer, chains, tol, seed)
    S = Set{Int}()                                   # vertices interior to an input edge
    foreach(c -> union!(S, c[2:end-1]), chains)
    setdiff!(S, outer)                               # outline corners stay corners
    nout = length(outer)
    whole = Set{NTuple{2,Int}}()                     # edges entered as they are
    pending = Vector{Int}[]                          # edges with an endpoint still to insert
    for (r, c) in enumerate(chains)
        r <= nout && continue                        # outline edges: the boundary curve
        c[1] in S || c[end] in S ? push!(pending, c) : push!(whole, minmax(c[1], c[end]))
    end
    tri = triangulate(pts; boundary_nodes = [outer; outer[1]],
        segments = isempty(whole) ? nothing : whole, skip_points = collect(S),
        rng = DT.Random.Xoshiro(seed))
    todo = copy(S)
    while !isempty(todo) || !isempty(pending)
        progress = false
        for v in collect(todo)                       # split the segment v lies on
            e = _segment_through(tri, pts, v, tol)
            e === nothing && continue
            DT.complete_split_edge_and_legalise!(tri, e[1], e[2], v)
            delete!(todo, v)
            progress = true
        end
        for c in collect(pending)                    # edges whose ends are now in place
            (c[1] in todo || c[end] in todo) && continue
            q = filter(v -> !(v in todo), c)
            for k in 1:length(q)-1
                DT.contains_segment(tri, q[k], q[k+1]) || DT.add_segment!(tri, q[k], q[k+1])
            end
            filter!(!=(c), pending)
            progress = true
        end
        if !progress                                 # cyclic overlaps of collinear edges
            v = pop!(todo)
            DT.add_point!(tri, v)
        end
    end
    return tri
end

# the constrained edge of `tri` that vertex v lies on, strictly between its ends
function _segment_through(tri, pts, v, tol)
    p = pts[v]
    for e in DT.each_segment(tri)
        i, j = DT.initial(e), DT.terminal(e)
        a, b = pts[i], pts[j]
        dx, dy = b[1] - a[1], b[2] - a[2]
        L = hypot(dx, dy)
        s = ((p[1] - a[1]) * dx + (p[2] - a[2]) * dy) / L
        abs((p[1] - a[1]) * dy - (p[2] - a[2]) * dx) / L <= tol && tol < s < L - tol && return (i, j)
    end
    return nothing
end

# expected triangle count of a graded mesh: ∫ dA / (√3/4·h²) over the material,
# by midpoint sampling on a 64 × 64 grid over the outline's bounding box
function _graded_count(h, outer, hs)
    xs, ys = first.(outer), last.(outer)
    x0, x1, y0, y1 = minimum(xs), maximum(xs), minimum(ys), maximum(ys)
    n = 64
    dx, dy = (x1 - x0) / n, (y1 - y0) / n
    total = 0.0
    for a in 1:n, b in 1:n
        x, y = x0 + (a - 0.5) * dx, y0 + (b - 0.5) * dy
        _inside(outer, x, y) && !any(H -> _inside(H, x, y), hs) || continue
        total += dx * dy / (sqrt(3) / 4 * Float64(h(x, y))^2)
    end
    return total
end

function _edge_nodes(mesh, a, b, tol)
    dx, dy = b[1] - a[1], b[2] - a[2]
    L = hypot(dx, dy)
    hits = Tuple{Float64,Int}[]
    for i in axes(mesh.points, 2)
        px, py = mesh.points[1, i] - a[1], mesh.points[2, i] - a[2]
        s = (px * dx + py * dy) / L
        -tol <= s <= L + tol && abs(px * dy - py * dx) / L <= tol && push!(hits, (s, i))
    end
    return last.(sort!(hits))
end

end # module
