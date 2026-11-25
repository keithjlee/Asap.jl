"""
    WallMesh{T, S}

A structured shell mesh of a planar rectangular wall with rectangular
openings, produced by [`mesh_wall`](@ref). The wall's own 2D coordinates are
`(s, z)`: `s` along its length (direction `along`), `z` up its height
(direction `up`); global position = `origin + s·along + z·up`.

# Fields
- `nodes::Vector{Node{T}}`: nodes used by the mesh (reused nodes included)
- `elements::Vector{ShellElement{T,4}}`: counter-clockwise about
  `along × up` (the wall's outward normal); local `e₁ = along`
- `s::Vector{T}`, `z::Vector{T}`: grid lines [length]
- `grid::Matrix{Union{Nothing,Node{T}}}`: node at grid point `(i, j)` (or
  `nothing` inside an opening)
- `cells::Matrix{Int}`: element index of grid cell `(i, j)` (0 = opening)
- `origin`, `along`, `up`: placement
"""
struct WallMesh{T,S<:AbstractShellSection{T}} <: AbstractShellMesh{T}
    nodes::Vector{Node{T}}
    elements::Vector{ShellElement{T,4,S}}
    s::Vector{T}
    z::Vector{T}
    grid::Matrix{Union{Nothing,Node{T}}}
    cells::Matrix{Int}
    origin::SVector{3,T}
    along::SVector{3,T}
    up::SVector{3,T}
end

# grid lines covering [0, L] through every breakpoint, intervals ≤ h
function _grid_lines(L, breaks, h)
    b = sort!(unique!(filter(x -> 0 <= x <= L, [0.0; L; breaks...])))
    lines = [b[1]]
    for k in 2:length(b)
        Δ = b[k] - b[k-1]
        Δ <= 1e-9 * L && continue
        n = max(1, ceil(Int, Δ / h - 1e-9))
        append!(lines, b[k-1] .+ Δ .* (1:n) ./ n)
    end
    return lines
end

"""
    mesh_wall(origin, along, up, length, height, section;
              openings = [], s_lines = [], z_lines = [], size = height / 4,
              reuse = Node[], fixity_bottom = :fixed, id = :wall) -> WallMesh

Mesh a planar rectangular wall with [`ShellElement`](@ref)s on a structured
grid.

- `origin`: global position of the wall's bottom-left corner
- `along`, `up`: unit in-plane directions (length, height); orthogonal
- `length`, `height` [length]
- `section`: the wall's [`AbstractShellSection`](@ref)
- `openings`: rectangles `(s₀, z₀, width, height)` in wall coordinates
  (doors, windows); cells inside are removed and opening edges become grid
  lines, so piers and spandrels are meshed exactly
- `s_lines`, `z_lines`: extra hard lines — storey levels, frame/column lines,
  load points — guaranteed to be grid lines
- `size`: target element size [length]
- `reuse`: existing nodes; a mesh node at the same position (to 1e-9
  relative) REUSES it, which is how walls connect to each other (L/C/T
  cores share their corner lines) and to frames (beam/column nodes on the
  wall's grid lines)
- `fixity_bottom`: support symbol for the bottom row (`:fixed`, `:pinned`,
  or `:free` to leave it unsupported); reused nodes keep their own fixity

The structured grid keeps every element rectangular — the regime where the
membrane is exact in pure bending (see `docs/shells/DERIVATION.md`).
"""
function mesh_wall(origin::AbstractVector, along::AbstractVector, up::AbstractVector,
    length::Real, height::Real, section::AbstractShellSection{T};
    openings = Tuple[], s_lines = Float64[], z_lines = Float64[], size::Real = height / 4,
    reuse::AbstractVector = Node{T}[], fixity_bottom::Symbol = :fixed, id::Symbol = :wall) where {T}
    o = SVector{3,T}(origin)
    a = normalize(SVector{3,T}(along))
    u = normalize(SVector{3,T}(up))
    abs(dot(a, u)) < 1e-9 || throw(ArgumentError("`along` and `up` must be orthogonal"))
    sb = Float64[s_lines...]; zb = Float64[z_lines...]
    for op in openings
        s0, z0, w, h = op
        (s0 >= 0 && z0 >= 0 && s0 + w <= length && z0 + h <= height && w > 0 && h > 0) ||
            throw(ArgumentError("opening $op lies outside the wall"))
        append!(sb, (s0, s0 + w)); append!(zb, (z0, z0 + h))
    end
    S = T.(_grid_lines(length, sb, size))
    Z = T.(_grid_lines(height, zb, size))
    ns, nz = Base.length(S), Base.length(Z)
    inopening(sc, zc) = any(op -> op[1] < sc < op[1] + op[3] && op[2] < zc < op[2] + op[4], openings)

    cellon = [!inopening((S[i] + S[i+1]) / 2, (Z[j] + Z[j+1]) / 2) for i in 1:ns-1, j in 1:nz-1]
    used = falses(ns, nz)
    for i in 1:ns-1, j in 1:nz-1
        cellon[i, j] && (used[i, j] = used[i+1, j] = used[i+1, j+1] = used[i, j+1] = true)
    end

    scale = max(length, height)
    keyof(x) = round.(Int, x ./ (1e-9 * scale))
    lookup = Dict{Any,Node{T}}(keyof(n.position) => n for n in reuse)
    grid = Matrix{Union{Nothing,Node{T}}}(nothing, ns, nz)
    nodes_ = Node{T}[]
    for i in 1:ns, j in 1:nz
        used[i, j] || continue
        x = o + S[i] * a + Z[j] * u
        n = get(lookup, keyof(x), nothing)
        if n === nothing
            fix = j == 1 && fixity_bottom !== :free ? fixity_bottom : :free
            n = Node(collect(x), fix, id)
            lookup[keyof(x)] = n
        end
        grid[i, j] = n
        push!(nodes_, n)
    end
    els = ShellElement{T,4,typeof(section)}[]
    cells = zeros(Int, ns - 1, nz - 1)
    for i in 1:ns-1, j in 1:nz-1
        cellon[i, j] || continue
        push!(els, ShellElement((grid[i, j], grid[i+1, j], grid[i+1, j+1], grid[i, j+1]), section, id))
        cells[i, j] = Base.length(els)
    end
    return WallMesh{T,typeof(section)}(nodes_, els, S, Z, grid, cells, o, a, u)
end

"""
    wall_cut(mesh::WallMesh, z; s_range = (-Inf, Inf)) -> (elements, cut_nodes)

The free body below height `z` (optionally restricted to `s_range` — a
pier between openings) and its nodes on the cut line, ready for
[`section_cut`](@ref). `z` must be a grid line (pass it in `z_lines`).
"""
function wall_cut(mesh::WallMesh{T}, z::Real; s_range = (-Inf, Inf)) where {T}
    j = findfirst(zz -> isapprox(zz, z; atol = 1e-9 * max(abs(z), one(T))), mesh.z)
    j === nothing && throw(ArgumentError("z = $z is not a grid line of this wall (add it to z_lines)"))
    lo, hi = s_range
    els = ShellElement{T}[]
    for i in axes(mesh.cells, 1), jj in 1:j-1
        c = mesh.cells[i, jj]
        c > 0 && lo <= (mesh.s[i] + mesh.s[i+1]) / 2 <= hi && push!(els, mesh.elements[c])
    end
    cut = Node{T}[]
    for i in eachindex(mesh.s)
        n = mesh.grid[i, j]
        n !== nothing && lo <= mesh.s[i] <= hi && any(el -> any(m -> m === n, el.nodes), els) && push!(cut, n)
    end
    return els, cut
end

"""
    edge_loads(el::ShellElement, edge::Int, w::AbstractVector; case = :LC1) -> Vector{NodeLoad}

Consistent nodal loads of a uniform LINE traction `w` [force/length]
(global vector) along edge `edge` (1…4; edge k joins nodes k → k+1) of a
4-node shell: forces `w·l/2` at both ends plus — because the drilling
membrane's edges deform quadratically with the end rotations — drilling
moments `∓(w·n)·l²/12` about the element normal at the start/end node,
where `n` is the edge's outward in-plane normal. Without the moments a
uniform edge traction is not represented consistently (e.g. a uniform
tension field would not be reproduced exactly).
"""
function edge_loads(el::ShellElement{T,4}, edge::Int, w::AbstractVector; case::Symbol = :LC1) where {T}
    1 <= edge <= 4 || throw(ArgumentError("edge must be 1…4"))
    nJ, nK = el.nodes[edge], el.nodes[mod1(edge + 1, 4)]
    e3 = local_frame(el)[3, :]
    d = nK.position - nJ.position
    l = norm(d)
    n = normalize(cross(d, e3))                 # outward for counter-clockwise nodes
    ws = SVector{3,T}(w)
    pn = dot(ws, n)
    F = ws * (l / 2)
    m = (pn * l^2 / 12) * SVector{3,T}(e3)
    return NodeLoad{T}[NodeForce(nJ, F; case = case), NodeForce(nK, F; case = case),
                       NodeMoment(nJ, -m; case = case), NodeMoment(nK, m; case = case)]
end
