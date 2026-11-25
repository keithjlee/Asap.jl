#=
Dependency-free reader for Gmsh meshes (MSH 4.1, ASCII): first-order triangles
and quadrilaterals become ShellElements, physical groups become named node and
element sets. Planar meshes keep in-plane coordinates (segment and region
queries); curved surfaces are read as spatial meshes.
=#

"""
    read_gmsh(path, section; planar = :auto, fixity = :free, id = :gmsh) -> ShellMesh

Read a Gmsh mesh file (format MSH 4.1, ASCII — Gmsh's default: `gmsh -2
model.geo -format msh41`) into a [`ShellMesh`](@ref) of [`ShellElement`](@ref)s.

- Every 3-node triangle and 4-node quadrilateral of a 2D entity becomes an
  element with `section` and tag `id`. Mesh first-order (`-order 1`); second-order
  elements are refused.
- **Physical groups** become `mesh.node_groups` (every group, any dimension:
  the nodes of a physical point, curve or surface — e.g. a support edge) and
  `mesh.element_groups` (physical surfaces), keyed by the group name (or
  `"dim:tag"` when unnamed). Read them with [`group_nodes`](@ref) and
  [`group_elements`](@ref).
- `planar = :auto`: a mesh whose nodes lie in one plane (within 1e-9 of its
  size) is planar: in-plane coordinates are kept (`mesh.points` is 2 × n; the
  plane's `e1` is the global X axis projected onto it, or Y if X is normal to
  it) and elements are ordered counter-clockwise about the plane's normal (the
  first element's normal). Otherwise the mesh is spatial (3 × n global
  coordinates) and elements keep Gmsh's orientation, which is consistent per
  surface. `planar = true` insists (error if not planar); `false` never projects.
- `fixity`: initial fixity of every node (`:free`); apply supports to groups
  with [`fixnode!`](@ref).

Nodes not used by any shell element (e.g. on 1D-only entities) are dropped.
"""
function read_gmsh(path::AbstractString, section::AbstractShellSection{T};
    planar = :auto, fixity::Symbol = :free, id::Symbol = :gmsh) where {T}
    planar in (:auto, true, false) || throw(ArgumentError("planar must be :auto, true or false"))
    lines = readlines(path)
    sec = _msh_sections(lines)
    haskey(sec, "MeshFormat") || throw(ArgumentError("$path: not a Gmsh mesh (no \$MeshFormat)"))
    fmt = split(lines[first(sec["MeshFormat"])])
    (startswith(fmt[1], "4.1") && fmt[2] == "0") || throw(ArgumentError(
        "$path: MSH $(fmt[1]) $(fmt[2] == "0" ? "ASCII" : "binary") is not supported; " *
        "export MSH 4.1 ASCII (gmsh … -format msh41, without -bin)"))
    # physical names and entity → physical tags
    names = Dict{Tuple{Int,Int},String}()
    if haskey(sec, "PhysicalNames")
        r = sec["PhysicalNames"]
        for l in lines[first(r)+1:last(r)]
            m = match(r"^\s*(\d+)\s+(\d+)\s+\"(.*)\"", l)
            m === nothing || (names[(parse(Int, m[1]), parse(Int, m[2]))] = m[3])
        end
    end
    phys = haskey(sec, "Entities") ? _msh_entities(lines, sec["Entities"]) : Dict{Tuple{Int,Int},Vector{Int}}()
    # nodes
    coords = Dict{Int,SVector{3,T}}()
    r = sec["Nodes"]; i = first(r)
    nblocks = parse(Int, split(lines[i])[1]); i += 1
    for _ in 1:nblocks
        _, _, par, nn = parse.(Int, split(lines[i])); i += 1
        tags = [parse(Int, lines[i+k]) for k in 0:nn-1]; i += nn
        for k in 1:nn
            v = parse.(T, split(lines[i+k-1]))
            coords[tags[k]] = SVector{3,T}(v[1], v[2], v[3])
        end
        i += nn
    end
    # elements: shell cells (dim 2) and node sets of every entity
    cells = Vector{Vector{Int}}()
    cell_entity = Tuple{Int,Int}[]
    entity_nodes = Dict{Tuple{Int,Int},Set{Int}}()
    r = sec["Elements"]; i = first(r)
    nblocks = parse(Int, split(lines[i])[1]); i += 1
    for _ in 1:nblocks
        dim, etag, etype, ne = parse.(Int, split(lines[i])); i += 1
        etype in (8, 9, 10, 16, 21, 23) && throw(ArgumentError(
            "$path: second-order elements (type $etype) are not supported; mesh with -order 1"))
        ns = get!(entity_nodes, (dim, etag), Set{Int}())
        for k in 0:ne-1
            v = parse.(Int, split(lines[i+k]))[2:end]
            union!(ns, v)
            if dim == 2
                etype in (2, 3) || throw(ArgumentError("$path: unsupported 2D element type $etype"))
                push!(cells, v); push!(cell_entity, (dim, etag))
            end
        end
        i += ne
    end
    isempty(cells) && throw(ArgumentError("$path: no triangles or quadrilaterals"))
    used = sort!(unique!(reduce(vcat, cells)))
    idx = Dict(t => k for (k, t) in enumerate(used))
    X = reduce(hcat, [coords[t] for t in used])
    rcells = [[idx[t] for t in c] for c in cells]
    # planarity
    Q, o, e1, e2 = _msh_plane(X, rcells, planar, path)
    m = _shell_mesh(T, Q, rcells, section, o, e1, e2, fixity, id)
    # physical groups → named sets
    ngroups = Dict{String,Vector{Int}}(); egroups = Dict{String,Vector{Int}}()
    for ((dim, etag), ptags) in phys, pt in ptags
        name = get(names, (dim, pt), "$dim:$pt")
        nodes_ = get(entity_nodes, (dim, etag), Set{Int}())
        append!(get!(ngroups, name, Int[]), [idx[t] for t in nodes_ if haskey(idx, t)])
        dim == 2 && append!(get!(egroups, name, Int[]), findall(==((dim, etag)), cell_entity))
    end
    foreach(v -> unique!(sort!(v)), values(ngroups)); foreach(v -> unique!(sort!(v)), values(egroups))
    return ShellMesh(m.nodes, m.elements, m.points, m.boundary, m.origin, m.e1, m.e2, ngroups, egroups)
end

# line ranges of every $Section … $EndSection block (contents only)
function _msh_sections(lines)
    out = Dict{String,UnitRange{Int}}()
    k = 1
    while k <= length(lines)
        l = strip(lines[k])
        if startswith(l, "\$") && !startswith(l, "\$End")
            name = l[2:end]
            stop = findnext(x -> strip(x) == "\$End" * name, lines, k + 1)
            stop === nothing && throw(ArgumentError("unterminated \$$name block"))
            out[name] = k+1:stop-1
            k = stop
        end
        k += 1
    end
    return out
end

# (dim, entity tag) → physical tags, from the $Entities block
function _msh_entities(lines, r)
    out = Dict{Tuple{Int,Int},Vector{Int}}()
    i = first(r)
    counts = parse.(Int, split(lines[i])); i += 1
    for (dim, n) in zip(0:3, counts), _ in 1:n
        v = split(lines[i]); i += 1
        tag = parse(Int, v[1])
        j = dim == 0 ? 5 : 8                      # point: x y z; others: bounding box
        np = parse(Int, v[j])
        np > 0 && (out[(dim, tag)] = parse.(Int, v[j+1:j+np]))
    end
    return out
end

# planar meshes → in-plane coordinates in a frame from the first element's normal
function _msh_plane(X::Matrix{T}, cells, planar, path) where {T}
    o = X[:, cells[1][1]]
    c = cells[1]
    n = cross(X[:, c[2]] - X[:, c[1]], X[:, c[3]] - X[:, c[1]])
    n = n / norm(n)
    scale = maximum(abs, X .- o)
    flat = maximum(k -> abs(dot(X[:, k] - o, n)), axes(X, 2)) <= 1e-9 * max(scale, one(T))
    if planar === false || (planar === :auto && !flat)
        return X, zeros(T, 3), T[1, 0, 0], T[0, 1, 0]
    end
    flat || throw(ArgumentError("$path: planar = true but the mesh is not planar"))
    ex = abs(n[1]) < 0.9 ? T[1, 0, 0] : T[0, 1, 0]
    e1 = ex - dot(ex, n) * n; e1 /= norm(e1)
    e2 = cross(n, e1)
    Q = reduce(hcat, [[dot(X[:, k] - o, e1), dot(X[:, k] - o, e2)] for k in axes(X, 2)])
    return Q, o, e1, e2
end
