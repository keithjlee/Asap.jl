#=
General planar shell meshes: an array constructor for meshes from any mesher
(triangles and quadrilaterals), boundary and segment queries, consistent edge
loads along boundary segments, and reaction totals over node groups — the
pieces a disturbed-region (strut-and-tie) analysis of a panel needs. The
polygon mesher `mesh_polygon` is provided by the DelaunayTriangulation.jl
extension (`AsapDelaunayExt`); this file only declares it.
=#

"""
    AbstractShellMesh{T}

Supertype of the shell mesh containers: [`WallMesh`](@ref) (structured walls
with rectangular openings) and [`ShellMesh`](@ref) (general planar meshes).
Every subtype carries `nodes::Vector{Node{T}}` and an element vector of
[`ShellElement`](@ref)s, ready for a `Model`.
"""
abstract type AbstractShellMesh{T} end

"""
    ShellMesh{T, E<:ShellElement} <: AbstractShellMesh{T}

A planar shell mesh in 3D space: nodes, elements and the boundary, with the
in-plane coordinates kept for geometric queries. Returned by
[`shell_mesh`](@ref) (from arrays of points and cells, e.g. another mesher's
output) and [`mesh_polygon`](@ref) (DelaunayTriangulation.jl extension).

# Fields
- `nodes::Vector{Node{T}}`: mesh nodes; node `i` sits at `points[:, i]`
- `elements::Vector{E}`: shell elements, counter-clockwise about `e₃`
- `points::Matrix{T}`: 2 × n in-plane coordinates of a planar mesh, or 3 × n
  global coordinates of a curved (spatial) one, e.g. from [`read_gmsh`](@ref)
  [length]
- `boundary::Vector{Vector{Vector{Int}}}`: boundary curves → edges → node
  indices in order. Curve 1 is the outline (counter-clockwise), curves 2… the
  holes (clockwise). For [`mesh_polygon`](@ref) each edge of an input polygon
  is its own entry, so `boundary[c][k]` lists every mesh node on input edge
  `k` of curve `c`, end to end; for [`shell_mesh`](@ref) each closed loop is
  a single entry.
- `origin`, `e1`, `e2`: the plane of a planar mesh: node `i` is at
  `origin + x·e1 + y·e2` (orthonormal `e1`, `e2`; element normals are
  `e3 = e1 × e2`); the global axes for a spatial mesh
- `node_groups`, `element_groups::Dict{String,Vector{Int}}`: named node and
  element index sets (Gmsh physical groups; empty otherwise)

Use [`boundary_nodes`](@ref), [`nodes_on_segment`](@ref), [`edge_loads`](@ref)
(segment form), [`elements_in`](@ref) (planar), [`group_nodes`](@ref),
[`group_elements`](@ref) and [`total_reaction`](@ref) to place supports and
loads and to read results by region.
"""
struct ShellMesh{T,E<:ShellElement} <: AbstractShellMesh{T}
    nodes::Vector{Node{T}}
    elements::Vector{E}
    points::Matrix{T}
    boundary::Vector{Vector{Vector{Int}}}
    origin::SVector{3,T}
    e1::SVector{3,T}
    e2::SVector{3,T}
    node_groups::Dict{String,Vector{Int}}
    element_groups::Dict{String,Vector{Int}}
end

ShellMesh(nodes, elements, points, boundary, origin, e1, e2) =
    ShellMesh(nodes, elements, points, boundary, origin, e1, e2,
        Dict{String,Vector{Int}}(), Dict{String,Vector{Int}}())

_planar(mesh::ShellMesh) = size(mesh.points, 1) == 2

"""
    group_nodes(mesh::ShellMesh, name) -> Vector{Node}

Nodes of the named group (a Gmsh physical group of any dimension: a point, the
nodes of a physical curve such as a support edge, or of a physical surface).
"""
group_nodes(mesh::ShellMesh, name::AbstractString) = mesh.nodes[_group(mesh.node_groups, name)]

"""
    group_elements(mesh::ShellMesh, name) -> Vector{ShellElement}

Elements of the named group (a Gmsh physical surface).
"""
group_elements(mesh::ShellMesh, name::AbstractString) = mesh.elements[_group(mesh.element_groups, name)]

_group(d, name) = haskey(d, name) ? d[name] :
    throw(ArgumentError("no group \"$name\"; available: $(join(sort!(collect(keys(d))), ", "))"))

function _plane_axes(::Type{T}, origin, e1, e2) where {T}
    o = SVector{3,T}(origin)
    a = normalize(SVector{3,T}(e1))
    b = SVector{3,T}(e2) - dot(SVector{3,T}(e2), a) * a
    norm(b) > sqrt(eps(T)) || throw(ArgumentError("e1 and e2 must not be parallel"))
    return o, a, normalize(b)
end

_points2(P::AbstractMatrix) = size(P, 1) == 2 ? P :
    throw(DimensionMismatch("points must be 2 × n (in-plane coordinates)"))
_points2(P::AbstractVector) = reduce(hcat, [collect(float.(p))[1:2] for p in P])

"""
    shell_mesh(points, cells, section; origin = [0, 0, 0], e1 = [1, 0, 0], e2 = [0, 1, 0],
               fixity = :free, id = :mesh, reuse = Node[]) -> ShellMesh

Build shell elements from a planar mesh given as arrays — the output of any
mesher (Gmsh, Triangle, a CAD export, a hand-made grid).

- `points`: 2 × n matrix (or vector of 2-tuples) of in-plane coordinates [length]
- `cells`: vector of 3- or 4-tuples (or vectors) of 1-based point indices;
  triangles and quadrilaterals may mix. Orientation is normalized: cells are
  reordered counter-clockwise in the plane, so element normals all point
  along `e₃ = e1 × e2`.
- `section`: the shell section of every element
- `origin`, `e1`, `e2`: the plane in global coordinates (`e2` is
  orthogonalized against `e1`)
- `fixity`: fixity of every NEW node (`:free` by default; place supports with
  [`fixnode!`](@ref) on [`nodes_on_segment`](@ref))
- `reuse`: existing nodes (frame joints, another wall's edge): a mesh node at the
  same position (to 1e-9 of the mesh size) IS that node, which is how a mesh
  connects to frames and other meshes. Reused nodes keep their fixity and appear
  in `mesh.nodes` too: `unique([frame_nodes; mesh.nodes])` is the tidy node list
  (a node listed twice is harmless, its spare slot is retired as inactive).

Points no cell uses are dropped. The boundary is recovered as closed loops of
edges used by a single cell.
"""
function shell_mesh(points, cells, section::AbstractShellSection{T};
    origin = zeros(3), e1 = [1.0, 0, 0], e2 = [0.0, 1, 0], fixity::Symbol = :free,
    id::Symbol = :mesh, reuse::AbstractVector = Node{T}[]) where {T}
    P = Matrix{T}(_pointsd(points))
    used = sort!(unique!(reduce(vcat, [collect(Int, c) for c in cells])))
    remap = Dict(old => new for (new, old) in enumerate(used))
    Q = P[:, used]
    cellsr = map(cells) do c
        ids = [remap[i] for i in c]
        length(ids) in (3, 4) || throw(ArgumentError("cells must have 3 or 4 nodes"))
        ids
    end
    return _shell_mesh(T, Q, cellsr, section, origin, e1, e2, fixity, id, reuse)
end

# planar (2 × n): positions in the given plane, cells made counter-clockwise
function _shell_mesh(::Type{T}, Q::Matrix{T}, cells, section, origin, e1, e2, fixity, id,
        reuse = Node{T}[]) where {T}
    if size(Q, 1) == 2
        o, a, b = _plane_axes(T, origin, e1, e2)
        X = [o + Q[1, i] * a + Q[2, i] * b for i in axes(Q, 2)]
        nodes = _mesh_nodes(X, reuse, fixity)
        oriented = [_signed_area2(Q, ids) < 0 ? reverse(ids) : ids for ids in cells]
        loops = _boundary_loops(oriented, l -> abs(_signed_area2(Q, l[1:end-1])), Q)
    else
        # spatial (3 × n): global positions; the cells keep the mesher's orientation
        # (a consistently oriented surface keeps consistent element normals)
        o, a, b = zero(SVector{3,T}), SVector{3,T}(1, 0, 0), SVector{3,T}(0, 1, 0)
        nodes = _mesh_nodes([SVector{3,T}(Q[1, i], Q[2, i], Q[3, i]) for i in axes(Q, 2)], reuse, fixity)
        oriented = cells
        loops = _boundary_loops(oriented, l -> sum(k -> norm(Q[:, l[k+1]] - Q[:, l[k]]), 1:length(l)-1))
    end
    els = [ShellElement(Tuple(nodes[ids]), section, id) for ids in oriented]
    return ShellMesh(nodes, identity.(els), Q, [[loop] for loop in loops], o, a, b)
end

# a node per position: an existing node from `reuse` where one sits there (to
# 1e-9 of the mesh size), else a new node with `fixity`
function _mesh_nodes(X::Vector{SVector{3,T}}, reuse, fixity) where {T}
    isempty(reuse) && return [Node(x, fixity) for x in X]
    scale = maximum(x -> maximum(abs, x), X)
    tol = 1e-9 * max(scale, one(T))
    return map(X) do x
        i = findfirst(n -> norm(n.position - x) <= tol, reuse)
        i === nothing ? Node(x, fixity) : reuse[i]
    end
end

_pointsd(P::AbstractMatrix) = size(P, 1) in (2, 3) ? P :
    throw(DimensionMismatch("points must be 2 × n (in-plane) or 3 × n (global)"))
_pointsd(P::AbstractVector) = reduce(hcat, [collect(float.(p)) for p in P])

_signed_area2(Q, ids) = sum(k -> Q[1, ids[k]] * Q[2, ids[mod1(k + 1, length(ids))]] -
                                 Q[1, ids[mod1(k + 1, length(ids))]] * Q[2, ids[k]], eachindex(ids))

# closed boundary loops (node indices, first == last) from oriented cells: an
# edge used by exactly one cell is a boundary edge, followed tail to head
function _boundary_loops(cells, size_of, Q = nothing)
    count = Dict{Tuple{Int,Int},Int}()
    for c in cells, k in eachindex(c)
        e = minmax(c[k], c[mod1(k + 1, length(c))])
        count[e] = get(count, e, 0) + 1
    end
    # boundary half-edges i → j (material on the left); a pinch vertex (a hole
    # touching the outline or another hole at a point) has several outgoing
    out = Dict{Int,Vector{Int}}()
    for c in cells, k in eachindex(c)
        i, j = c[k], c[mod1(k + 1, length(c))]
        count[minmax(i, j)] == 1 && push!(get!(out, i, Int[]), j)
    end
    loops = Vector{Vector{Int}}()
    while !isempty(out)
        start = minimum(keys(out))
        loop = [start]
        prev, i = 0, start
        while true
            j = _next_boundary(out, prev, i, Q)
            push!(loop, j)
            prev, i = i, j
            j == start && !haskey(out, start) && break
            j == start && length(loop) > 1 && break
        end
        push!(loops, loop)
    end
    # outline first: the largest loop (planar: enclosed area, counter-clockwise;
    # holes come out clockwise because cells are counter-clockwise)
    sort!(loops; by = l -> -size_of(l))
    return loops
end

# take (and remove) the outgoing boundary edge at i: the only one, or at a pinch
# the first one clockwise from the incoming direction i → prev (keeps the material
# on the left and splits a pinched boundary into simple loops)
function _next_boundary(out, prev, i, Q)
    cand = out[i]
    k = 1
    if length(cand) > 1 && Q !== nothing && size(Q, 1) == 2 && prev != 0
        back = atan(Q[2, prev] - Q[2, i], Q[1, prev] - Q[1, i])
        cw(j) = mod(back - atan(Q[2, j] - Q[2, i], Q[1, j] - Q[1, i]), 2π)
        k = argmin(j -> cw(cand[j]) == 0 ? 2π : cw(cand[j]), eachindex(cand))
    end
    j = cand[k]
    deleteat!(cand, k)
    isempty(cand) && delete!(out, i)
    return j
end

"""
    boundary_nodes(mesh::ShellMesh, curve = 1, edge = nothing) -> Vector{Node}

Nodes on boundary `curve` (1 = outline, 2… = holes), in order: the whole
closed curve (`edge = nothing`, first node not repeated), or only input edge
`edge` of that curve (end nodes included) for a [`mesh_polygon`](@ref) mesh.
"""
function boundary_nodes(mesh::ShellMesh, curve::Int = 1, edge::Union{Nothing,Int} = nothing)
    edges = mesh.boundary[curve]
    if edge === nothing
        ids = Int[]
        for e in edges
            append!(ids, isempty(ids) ? e : e[2:end])
        end
        length(ids) > 1 && ids[end] == ids[1] && pop!(ids)
        return mesh.nodes[ids]
    end
    return mesh.nodes[edges[edge]]
end

"""
    nodes_on_segment(mesh::ShellMesh, a, b; tol = 1e-8·scale) -> Vector{Node}

Mesh nodes lying on the in-plane segment from `a` to `b` (2-vectors in the
mesh's plane coordinates), ordered from `a` to `b` — e.g. the nodes of a
bearing or of a loaded edge. Make sure the mesh has vertices at `a` and `b`
(polygon corners or extra points of [`mesh_polygon`](@ref)).
"""
function nodes_on_segment(mesh::ShellMesh{T}, a, b; tol = nothing) where {T}
    d_ = size(mesh.points, 1)
    length(a) == d_ == length(b) || throw(DimensionMismatch(
        "segment end points must be $(d_)-vectors (the mesh's $(d_ == 2 ? "in-plane" : "global") coordinates)"))
    A = collect(T, a)
    d = collect(T, b) - A
    L = norm(d)
    L > 0 || throw(ArgumentError("segment end points coincide"))
    scale = maximum(abs, mesh.points)
    τ = tol === nothing ? 1e-8 * max(scale, L) : tol
    hits = Tuple{T,Int}[]
    P = mesh.points
    for i in axes(P, 2)
        # distance along the segment, and off it (no temporaries)
        t = sum(k -> (P[k, i] - A[k]) * d[k], 1:d_) / L
        -τ <= t <= L + τ || continue
        off2 = sum(k -> (P[k, i] - A[k] - t * d[k] / L)^2, 1:d_)
        off2 <= τ^2 && push!(hits, (t, i))
    end
    sort!(hits)
    return mesh.nodes[last.(hits)]
end

"""
    edge_loads(mesh::ShellMesh, a, b, w; case = :LC1) -> Vector{NodeLoad}

Consistent nodal loads of a uniform line traction `w` [force/length] (global
3-vector) along the boundary segment `a → b` of a mesh (in-plane
coordinates): the element-edge loads ([`edge_loads`](@ref)) of every
boundary element edge whose end nodes both lie on the segment, so forces and
the drilling moments a drilling membrane needs are both included.
"""
function edge_loads(mesh::ShellMesh{T}, a, b, w::AbstractVector; case::Symbol = :LC1) where {T}
    on = Set(nodes_on_segment(mesh, a, b))
    loads = NodeLoad{T}[]
    for el in mesh.elements, k in 1:length(el.nodes)
        nJ, nK = el.nodes[k], el.nodes[mod1(k + 1, length(el.nodes))]
        (nJ in on && nK in on) || continue
        append!(loads, edge_loads(el, k, w; case = case))
    end
    isempty(loads) && throw(ArgumentError("no element edge lies on the segment $(a) → $(b)"))
    return loads
end

"""
    edge_loads(el::ShellElement{T,3}, edge, w; case = :LC1) -> Vector{NodeLoad}

Triangle version (edges 1…3, edge k joins nodes k → k+1): forces `w·l/2` at
both ends and drilling moments `∓α_b·(w·n)·l²/12` with α_b = 3/2, the
ANDES-OPT force-lumping parameter (Felippa 2003), where `n` is the edge's
outward in-plane normal. With these moments a uniform traction reproduces a
uniform stress field exactly (patch test).
"""
function edge_loads(el::ShellElement{T,3}, edge::Int, w::AbstractVector; case::Symbol = :LC1) where {T}
    1 <= edge <= 3 || throw(ArgumentError("edge must be 1…3"))
    nJ, nK = el.nodes[edge], el.nodes[mod1(edge + 1, 3)]
    e3 = local_frame(el)[3, :]
    d = nK.position - nJ.position
    l = norm(d)
    n = normalize(cross(d, e3))
    ws = SVector{3,T}(w)
    F = ws * (l / 2)
    m = (_ANDES_ALPHA_B * dot(ws, n) * l^2 / 12) * SVector{3,T}(e3)
    return NodeLoad{T}[NodeForce(nJ, F; case = case), NodeForce(nK, F; case = case),
                       NodeMoment(nJ, -m; case = case), NodeMoment(nK, m; case = case)]
end

"""
    elements_in(mesh::ShellMesh, polygon) -> Vector{ShellElement}

Elements whose centroid lies inside `polygon` (vector of in-plane 2-points or
2 × n matrix) — e.g. the patch under a bearing plate for an
[`AreaLoad`](@ref), or a region to read stresses from.
"""
function elements_in(mesh::ShellMesh{T}, polygon) where {T}
    _planar(mesh) || throw(ArgumentError("elements_in needs a planar mesh (in-plane coordinates)"))
    P = _points2(polygon)
    idx = Dict(n => i for (i, n) in enumerate(mesh.nodes))
    return filter(mesh.elements) do el
        c = sum(n -> mesh.points[:, idx[n]], el.nodes) / length(el.nodes)
        _inside(P, c[1], c[2])
    end
end

# even–odd ray casting
function _inside(P, x, y)
    inside = false
    n = size(P, 2)
    for k in 1:n
        xi, yi = P[1, k], P[2, k]
        xj, yj = P[1, mod1(k + 1, n)], P[2, mod1(k + 1, n)]
        ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) && (inside = !inside)
    end
    return inside
end

"""
    total_reaction(results, nodes; about = centroid) -> (F, M)

Sum of the support reactions over a group of nodes — e.g. one bearing of a
panel: the resultant force `F` and its moment `M` about `about` (default: the
centroid of the nodes), global axes. Works with any `LinearResults` (a case
or a combination).
"""
function total_reaction(res::LinearResults{T}, nodes::AbstractVector; about = nothing) where {T}
    c = about === nothing ? sum(n -> n.position, nodes) / length(nodes) : SVector{3,T}(about)
    F = zero(SVector{3,T}); M = zero(SVector{3,T})
    for n in nodes
        r = reaction(res, n)
        f = SVector{3,T}(r[1], r[2], r[3])
        F += f
        M += SVector{3,T}(r[4], r[5], r[6]) + cross(n.position - c, f)
    end
    return F, M
end

"""
    circle_polygon(center, radius; n = 32) -> Vector{NTuple{2}}

A regular `n`-gon inscribed in a circle (counter-clockwise), for circular
openings in [`mesh_polygon`](@ref). Choose `n` from a chord tolerance δ:
`n ≥ π / acos(1 − δ/radius)`.
"""
circle_polygon(center, radius::Real; n::Int = 32) =
    [(center[1] + radius * cospi(2k / n), center[2] + radius * sinpi(2k / n)) for k in 0:n-1]

"""
    mesh_polygon(outline, section; holes = [], h = nothing, max_area = nothing,
                 min_angle = 25, points = [], origin = [0, 0, 0], e1 = [1, 0, 0],
                 e2 = [0, 1, 0], fixity = :free, id = :panel, seed = 1, reuse = Node[]) -> ShellMesh

Triangulate a planar polygon with polygonal holes into triangular
[`ShellElement`](@ref)s: constrained Delaunay with quality refinement
(minimum angle `min_angle` [degrees]; maximum triangle area `max_area`, or a
target edge length `h`). **Requires DelaunayTriangulation.jl**
(`using DelaunayTriangulation` activates the extension).

- `outline`: polygon vertices in plane coordinates (vector of 2-points or
  2 × n), any orientation; include as vertices every point where a support or
  a load starts or ends, so meshes nodes land there
- `holes`: polygons inside the outline (any orientation); circles via
  [`circle_polygon`](@ref). Holes may TOUCH the outline or each other — at a
  vertex, with a vertex on an edge, or along a shared edge (a notch drawn as a
  hole, adjacent openings drawn as two) — but must not cross the outline or
  overlap each other (refused with an `ArgumentError`)
- `points`: extra interior points that must become nodes (point loads)
- `origin`, `e1`, `e2`, `fixity`, `id`: as in [`shell_mesh`](@ref)
- `seed`: the triangulator inserts points in random order; the seeded generator
  makes the mesh reproducible (same input and seed → same mesh)
- `reuse`: existing nodes (frame joints, beam and column nodes, another wall's
  edge nodes). Each one lying in the panel's plane, inside the material or on its
  boundary, becomes a mesh vertex and the mesh uses THAT node object, which
  connects the panel to them; nodes off the panel are ignored. Reused nodes keep
  their fixity and appear in `mesh.nodes` too: `unique([frame_nodes; mesh.nodes])`
  is the tidy node list (a node listed twice is harmless).

Each polygon edge is tracked through refinement: `mesh.boundary[c][k]` holds
every node on edge `k` of curve `c` (1 = outline, 1 + j = hole j), and
[`boundary_nodes`](@ref)`(mesh, c, k)` returns them.
"""
function mesh_polygon end
