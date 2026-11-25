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
cannot express. After refinement the triangles outside the material (outside
the outline or inside a hole) are dropped, and the nodes on every input edge
are recovered geometrically.
"""
module AsapDelaunayExt

using Asap
using Asap: ShellMesh, AbstractShellSection
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
    holes = [], h = nothing, max_area = nothing, min_angle::Real = 25.0, points = [],
    origin = zeros(3), e1 = [1.0, 0, 0], e2 = [0.0, 1, 0], fixity::Symbol = :free,
    id::Symbol = :panel, seed::Integer = 1) where {T}
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
    foreach(p -> vid((Float64(p[1]), Float64(p[2]))), points)
    # every input edge, split where another vertex lies on it (T-junctions); an
    # edge shared by two polygons (a hole along the outline or another hole) once
    edges = Set{NTuple{2,Int}}()
    for poly in ids, k in eachindex(poly)
        i, j = poly[k], poly[mod1(k + 1, length(poly))]
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
        for q in 1:length(chain)-1
            push!(edges, minmax(chain[q], chain[q+1]))
        end
    end
    # crossing segments: a hole crosses the outline or two holes overlap
    E = collect(edges)
    for x in eachindex(E), y in x+1:length(E)
        (a, b), (c, d) = E[x], E[y]
        _crosses(pts[a], pts[b], pts[c], pts[d], tol * scale) && throw(ArgumentError(
            "mesh_polygon: polygon edges cross near $(pts[a]) → $(pts[b]): holes must lie inside " *
            "the outline and must not overlap each other (touching is allowed)"))
    end
    # triangulate the segment graph (seeded: the same input gives the same mesh)
    tri = triangulate(pts; segments = edges, rng = DT.Random.Xoshiro(seed))
    amax = max_area !== nothing ? Float64(max_area) :
           h !== nothing ? sqrt(3) / 4 * Float64(h)^2 : nothing
    if amax === nothing
        refine!(tri; min_angle = Float64(min_angle), rng = DT.Random.Xoshiro(seed))
    else
        refine!(tri; min_angle = Float64(min_angle), max_area = amax, rng = DT.Random.Xoshiro(seed))
    end
    # keep the material: inside the outline, outside every hole (centroid test)
    P = DT.get_points(tri)
    cells = Vector{Int}[]
    for t in DT.each_solid_triangle(tri)
        i, j, k = DT.triangle_vertices(t)
        cx, cy = (P[i][1] + P[j][1] + P[k][1]) / 3, (P[i][2] + P[j][2] + P[k][2]) / 3
        _inside(outer, cx, cy) && !any(H -> _inside(H, cx, cy), hs) && push!(cells, [i, j, k])
    end
    mesh = Asap.shell_mesh([P[i] for i in eachindex(P)], cells, section;
        origin = origin, e1 = e1, e2 = e2, fixity = fixity, id = id)
    # the mesh nodes on every input edge, in order along it
    boundary = [[_edge_nodes(mesh, poly[k], poly[mod1(k + 1, length(poly))], tol) for k in eachindex(poly)]
                for poly in polys]
    return ShellMesh(mesh.nodes, mesh.elements, mesh.points, boundary, mesh.origin, mesh.e1, mesh.e2)
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
