"""
    AsapDelaunayExt

`mesh_polygon`: planar polygons with polygonal holes → triangular shell
meshes, through DelaunayTriangulation.jl (pure Julia: constrained Delaunay
triangulation with Ruppert-type quality refinement). Loaded when
DelaunayTriangulation is in the environment; core Asap stays dependency-free.

Every edge of the input polygons is passed as its own boundary section, so
refinement keeps the list of nodes on each input edge — exactly the node sets
that supports and edge loads are applied to.
"""
module AsapDelaunayExt

using Asap
using Asap: ShellMesh, AbstractShellSection, SVector
using DelaunayTriangulation
const DT = DelaunayTriangulation

# 2D points of a polygon as a vector of Float64 tuples, no closing repeat
function _poly(P)
    pts = P isa AbstractMatrix ? [(Float64(P[1, k]), Float64(P[2, k])) for k in axes(P, 2)] :
          [(Float64(p[1]), Float64(p[2])) for p in P]
    length(pts) > 1 && pts[end] == pts[1] && pop!(pts)
    length(pts) >= 3 || throw(ArgumentError("a polygon needs at least 3 vertices"))
    return pts
end

_area2(p) = sum(k -> p[k][1] * p[mod1(k + 1, length(p))][2] - p[mod1(k + 1, length(p))][1] * p[k][2], eachindex(p))

function Asap.mesh_polygon(outline, section::AbstractShellSection{T};
    holes = [], h = nothing, max_area = nothing, min_angle::Real = 25.0, points = [],
    origin = zeros(3), e1 = [1.0, 0, 0], e2 = [0.0, 1, 0], fixity::Symbol = :free,
    id::Symbol = :panel) where {T}
    outer = _poly(outline)
    _area2(outer) < 0 && reverse!(outer)                    # outline counter-clockwise
    hs = map(holes) do H
        p = _poly(H)
        _area2(p) > 0 ? reverse(p) : p                       # holes clockwise
    end
    pts = Tuple{Float64,Float64}[]
    curves = Vector{Vector{Vector{Int}}}()
    for poly in (outer, hs...)
        start = length(pts)
        append!(pts, poly)
        n = length(poly)
        push!(curves, [[start + k, start + mod1(k + 1, n)] for k in 1:n])
    end
    extra = [(Float64(p[1]), Float64(p[2])) for p in points]
    append!(pts, extra)
    tri = triangulate(pts; boundary_nodes = curves)
    amax = max_area !== nothing ? Float64(max_area) :
           h !== nothing ? sqrt(3) / 4 * Float64(h)^2 : nothing
    if amax === nothing
        refine!(tri; min_angle = Float64(min_angle))
    else
        refine!(tri; min_angle = Float64(min_angle), max_area = amax)
    end
    P = DT.get_points(tri)
    cells = [collect(DT.triangle_vertices(t)) for t in DT.each_solid_triangle(tri)]
    # used points only, renumbered in first-use order of the cells
    mesh = Asap.shell_mesh([P[i] for i in eachindex(P)], cells, section;
        origin = origin, e1 = e1, e2 = e2, fixity = fixity, id = id)
    # map DT point indices to mesh node indices through the coordinates
    lookup = Dict{Tuple{Float64,Float64},Int}()
    for i in axes(mesh.points, 2)
        lookup[(Float64(mesh.points[1, i]), Float64(mesh.points[2, i]))] = i
    end
    boundary = [[[lookup[(P[i][1], P[i][2])] for i in DT.get_boundary_nodes(tri, c, k)]
                 for k in eachindex(curves[c])] for c in eachindex(curves)]
    return ShellMesh(mesh.nodes, mesh.elements, mesh.points, boundary, mesh.origin, mesh.e1, mesh.e2)
end

end # module
