"""
    Node{T<:Real}

A structural node: a point in 3D space carrying six degrees of freedom —
three translations (along global X, Y, Z) and three rotations (about global
X, Y, Z).

Every node exposes all six DOF slots regardless of what connects to it.
Which slots actually participate in a solve is decided later by *activity
analysis*: a DOF touched by no element stiffness (e.g. the rotations of a
node connected only to truss elements) is excluded automatically. This is
what lets truss and frame elements mix freely in one model — there is no
separate truss node type.

A `Node` is pure definition data. Analysis results (displacements,
reactions) live in the results object returned by the solve, accessed as
`displacement(results, node)` / `reaction(results, node)` — never on the
node itself.

# Fields
- `position::SVector{3,T}`: coordinates in the global system [length]
- `fixity::SVector{6,Bool}`: support condition per DOF, ordered
  (Tx, Ty, Tz, Rx, Ry, Rz). `true` = **free** to move, `false` = fixed by a
  support (legacy convention preserved).
- `id::Symbol`: user tag for group queries (e.g. `:support`, `:topchord`)
- `index::Int`: position of this node in its model's node vector; assigned
  by `process!`, `0` until then. Internal bookkeeping — not user-set.

# Constructors
    Node(position, fixity::Symbol, id = :node)         # common support types
    Node(position, fixity::AbstractVector{Bool}, id = :node)  # explicit DOFs
    Node{T}(position, fixity, id)                       # explicit scalar type

`position` may be any 3-element vector; it is converted to an `SVector{3,T}`.

Available fixity symbols (see [`FIXITIES`](@ref)):

| symbol | meaning |
|---|---|
| `:free` | all six DOFs free |
| `:fixed` | fully clamped (all six fixed) |
| `:pinned` | translations fixed, rotations free |
| `:xfixed`/`:yfixed`/`:zfixed` | one translation fixed, all else free |
| `:xfree`/`:yfree`/`:zfree` | one translation free, all else fixed |

# Examples
```julia-repl
julia> Node([0.0, 0.0, 3.0], :pinned)

julia> Node([2.5, 0.0, 3.0], :free, :midspan)

julia> Node([5.0, 0.0, 0.0], [false, false, false, true, true, true])  # ≡ :pinned
```

See also [`fixnode!`](@ref), [`planarize!`](@ref).
"""
mutable struct Node{T<:Real}
    position::SVector{3,T}
    fixity::SVector{6,Bool}
    id::Symbol
    index::Int

    function Node{T}(position::AbstractVector{<:Real}, fixity::AbstractVector{Bool},
        id::Symbol=:node) where {T<:Real}
        @assert length(position) == 3 "node position must be a 3-vector (global X, Y, Z)"
        @assert length(fixity) == 6 "fixity must have 6 entries (Tx, Ty, Tz, Rx, Ry, Rz)"
        return new{T}(SVector{3,T}(position), SVector{6,Bool}(fixity), id, 0)
    end
end

Node(position::AbstractVector{T}, fixity::AbstractVector{Bool}, id::Symbol=:node) where {T<:Real} =
    Node{float(T)}(position, fixity, id)

function Node(position::AbstractVector{<:Real}, fixity::Symbol, id::Symbol=:node)
    @assert haskey(FIXITIES, fixity) "unknown fixity :$fixity — choose from $(sort!(collect(keys(FIXITIES))))"
    return Node(position, FIXITIES[fixity], id)
end

Node{T}(position::AbstractVector{<:Real}, fixity::Symbol, id::Symbol=:node) where {T<:Real} =
    Node{T}(position, FIXITIES[fixity], id)

"""
    FIXITIES

Map from support-condition symbols to 6-entry DOF fixity vectors, ordered
(Tx, Ty, Tz, Rx, Ry, Rz) with `true` = free. Matches the legacy `fixDict`
exactly, so existing model-building code keeps its meaning.
"""
const FIXITIES = Dict{Symbol,SVector{6,Bool}}(
    :free => SVector(true, true, true, true, true, true),
    :fixed => SVector(false, false, false, false, false, false),
    :pinned => SVector(false, false, false, true, true, true),
    :xfixed => SVector(false, true, true, true, true, true),
    :yfixed => SVector(true, false, true, true, true, true),
    :zfixed => SVector(true, true, false, true, true, true),
    :xfree => SVector(true, false, false, false, false, false),
    :yfree => SVector(false, true, false, false, false, false),
    :zfree => SVector(false, false, true, false, false, false),
)

"""
    fixnode!(node::Node, fixity::Symbol)

Replace a node's support condition with one of the standard types in
[`FIXITIES`](@ref) (e.g. `:pinned`, `:fixed`, `:zfixed`).
"""
function fixnode!(node::Node, fixity::Symbol)
    @assert haskey(FIXITIES, fixity) "unknown fixity :$fixity — choose from $(sort!(collect(keys(FIXITIES))))"
    node.fixity = FIXITIES[fixity]
    return node
end

# DOF indices fixed by planarization: for each working plane, the
# out-of-plane translation and the two in-plane-axis rotations.
const PLANES = Dict{Symbol,SVector{3,Int}}(
    :XY => SVector(3, 4, 5),
    :YZ => SVector(1, 5, 6),
    :ZX => SVector(2, 4, 6),
)

"""
    planarize!(nodes, plane = :XY)

Constrain nodes to behave as a 2D structure in the given global plane by
fixing the out-of-plane DOFs: the translation normal to the plane and the
two rotations whose axes lie in the plane (e.g. for `:XY`: Tz, Rx, Ry).

Operates on a vector of nodes or a single node.
"""
function planarize!(node::Node, plane::Symbol=:XY)
    @assert haskey(PLANES, plane) "unknown plane :$plane — choose from :XY, :YZ, :ZX"
    fixity = MVector(node.fixity)
    fixity[PLANES[plane]] .= false
    node.fixity = SVector(fixity)
    return node
end

function planarize!(nodes::Vector{<:Node}, plane::Symbol=:XY)
    for node in nodes
        planarize!(node, plane)
    end
    return nodes
end
