"""
    ShellElement{T, N, S<:AbstractShellSection{T}} <: AbstractElement{T}

A flat shell (surface) element with `N` corner nodes — membrane action with a
true drilling rotation plus plate bending — for walls, slabs, diaphragms and
general shell surfaces. Mixes freely with frames and trusses in one `Model`;
shared nodes connect them with full 6-DOF continuity (the drilling DOF lets
an in-plane beam transfer moment into a wall).

Two shapes, mixing freely in one mesh:

- `N = 4`, quadrilateral: the Ibrahimbegović–Taylor–Wilson drilling membrane
  combined with the MITC4 plate.
- `N = 3`, triangle: Felippa's ANDES-OPT drilling membrane combined with the
  MITC3+ plate. Both plates are Reissner–Mindlin (thick and thin), so mixed
  quad/triangle meshes use one plate theory.

The node count is a type parameter so elements of different shapes assemble
in separate type-stable groups.

# Fields
- `nodes::NTuple{N,Node{T}}`: corner nodes, ordered **counter-clockwise**
  about the element's outward normal (right-hand rule); the normal sets the
  sign of the local z-axis, of out-of-plane loads, and of moment resultants
- `section::S`: the shell section ([`ShellSection`](@ref) or
  [`ShellRigiditySection`](@ref))
- `id::Symbol`: user tag
- `index::Int`: position in `model.elements` (assigned by `process!`)
- `drilling::NTuple{N,Float64}`: drilling force-lumping parameter of each
  edge (edge k joins node k → k+1) [-] (assigned by `process!`). Quadrilateral:
  1 (the ITW membrane's lumping, fixed). Triangle: α_b = 3/2 (ANDES-OPT), but 1
  on an edge it shares with a quadrilateral, so the two membranes lump a
  constant stress state alike and mixed meshes pass the patch test.

# Local axes
`e₃` = normalized cross product of the diagonals (quad) or of edges 1 and 3
(triangle); `e₁` = edge node 1 → node 2 projected onto the element plane;
`e₂ = e₃ × e₁` (see [`shell_frame`](@ref)).
Stress resultants are reported in these axes — order the nodes so edge 1 runs
along the direction you want as "1" (the wall mesher does this
consistently).

# Constructors
    ShellElement(nodes, section, id = :shell)   # nodes: tuple or vector of 3 or 4 Nodes

# Examples
```julia-repl
julia> wall = ShellSection(Material(4.0e6, 0.15, 0.2), 8 / 12);

julia> el = ShellElement((n1, n2, n3, n4), wall, :wall)
```
"""
mutable struct ShellElement{T,N,S<:AbstractShellSection{T}} <: AbstractElement{T}
    nodes::NTuple{N,Node{T}}
    section::S
    id::Symbol
    index::Int
    drilling::NTuple{N,Float64}

    function ShellElement(nodes::NTuple{N,Node{T}}, section::S,
        id::Symbol = :shell) where {T,N,S<:AbstractShellSection{T}}
        N in (3, 4) || throw(ArgumentError(
            "ShellElement with $N nodes is not implemented (triangles: N = 3, quadrilaterals: N = 4)"))
        allunique(nodes) || throw(ArgumentError("shell element nodes must be distinct"))
        return new{T,N,S}(nodes, section, id, 0, N == 3 ? _ANDES_DRILL : ntuple(_ -> 1.0, Val(N)))
    end
end

ShellElement(nodes::AbstractVector{<:Node}, section::AbstractShellSection, id::Symbol = :shell) =
    ShellElement(Tuple(nodes), section, id)

nodes(el::ShellElement) = el.nodes

"""
    dof_signature(el::ShellElement) -> NTuple{6N, Bool}

All six DOFs of every node: the membrane stiffens the in-plane translations
AND the drilling rotation (a true drilling formulation, not a fictitious
spring), and the plate stiffens the out-of-plane translation and both
bending rotations. Every slot is active in any orientation.
"""
dof_signature(::ShellElement{T,N}) where {T,N} = ntuple(_ -> true, 6 * N)

"""
    stiffness(el::ShellElement{T,4}, x1, x2, x3, x4) -> SMatrix{24,24}

Global element stiffness at the given node positions (see
[`shell_stiffness`](@ref)).
"""
stiffness(el::ShellElement{T,4}, x1, x2, x3, x4) where {T} = shell_stiffness(el.section, x1, x2, x3, x4)
stiffness(el::ShellElement{T,3}, x1, x2, x3) where {T} = _shell_stiffness_d(el.section, (x1, x2, x3), el.drilling)

"""
    local_frame(el::ShellElement) -> SMatrix{3,3}

Rows are the element's local axes `e₁, e₂, e₃` in global coordinates.
"""
local_frame(el::ShellElement) = shell_frame(map(n -> n.position, el.nodes)...)[1]

# Per-edge drilling force lumping of the triangles (called by `process!`, after
# node indexing). A constant stress state puts drilling moments ∓α·t_n·ℓ²/12 at
# the ends of every element edge: the ITW quad lumps with α = 1, the ANDES-OPT
# triangle with α_b = 3/2. Across an edge shared by two elements the forces
# cancel, but the moments cancel only if both use the same α, so a triangle
# takes α = 1 on every edge a quadrilateral shares (also at folds and
# non-manifold junctions: matching α removes the in-plane part of the mismatch
# and is neutral otherwise). All other triangle edges keep α_b = 3/2.
function _drilling_lumping!(elements)
    quad_edges = Set{NTuple{2,Int}}()
    for el in elements
        el isa ShellElement{<:Any,4} || continue
        for k in 1:4
            push!(quad_edges, minmax(el.nodes[k].index, el.nodes[mod1(k + 1, 4)].index))
        end
    end
    for el in elements
        el isa ShellElement{<:Any,3} || continue
        el.drilling = ntuple(Val(3)) do k
            minmax(el.nodes[k].index, el.nodes[mod1(k + 1, 3)].index) in quad_edges ? 1.0 : _ANDES_DRILL[k]
        end
    end
    return elements
end

"""
    area(el::ShellElement) -> T

Element area [length²] (of the projection onto the mean plane), integrated
exactly with the bilinear map.
"""
function area(el::ShellElement{T,4}) where {T}
    _, xy, _ = shell_frame(map(n -> n.position, el.nodes)...)
    return _q4_area(xy)
end
function area(el::ShellElement{T,3}) where {T}
    _, xy, _ = shell_frame(map(n -> n.position, el.nodes)...)
    return _tri_area(xy)
end

# exact bilinear area: shoelace formula on the projected corners
_q4_area(xy) = sum(k -> xy[k][1] * xy[mod1(k + 1, 4)][2] - xy[mod1(k + 1, 4)][1] * xy[k][2], 1:4) / 2

_element_volume(el::ShellElement{T}) where {T} =
    el.section isa ShellSection ? thickness(el.section) * area(el) : zero(T)

"""
Maximum warping ratio (out-of-plane offset / √area) beyond which a 4-node flat
shell is rejected at `process!`: the flat-element projection error grows with
warping, and a strongly warped "quad" is usually a node-ordering mistake.
Mildly warped elements (≤ this ratio) are accepted (projected onto the mean
plane); a warping correction is planned (S7).
"""
const SHELL_WARP_LIMIT = 0.1

function validate_rigidities(el::ShellElement{T,4}) where {T}
    R, xy, warp = shell_frame(map(n -> n.position, el.nodes)...)
    warp <= SHELL_WARP_LIMIT || throw(ArgumentError(
        "shell element $(el.index) (:$(el.id)) is warped beyond the limit " *
        "(ratio $(round(warp; sigdigits = 3)) > $(SHELL_WARP_LIMIT)); check the node order " *
        "or split it into flatter elements"))
    # positive Jacobian at every corner ⇔ convex quad with counter-clockwise nodes
    for a in 1:4
        _, Nξ, Nη = _q4_N(_Q4_ξ[a], _Q4_η[a])
        detJ = _q4_jacobian(xy, Nξ, Nη)[5]
        detJ > 0 || throw(ArgumentError(
            "shell element $(el.index) (:$(el.id)) is not a convex quadrilateral with " *
            "consistently ordered nodes (non-positive Jacobian at node $a)"))
    end
    return nothing
end

# triangles are always flat; reject (near-)degenerate ones, whose frame and
# kernels are undefined: collinear nodes or a sliver below 1e-10 of the
# squared longest edge
function validate_rigidities(el::ShellElement{T,3}) where {T}
    x = map(n -> _sv3(n.position), el.nodes)
    a2 = norm(cross(x[2] - x[1], x[3] - x[1]))
    lmax = maximum(((1, 2), (2, 3), (3, 1))) do (i, j)
        norm(x[j] - x[i])
    end
    a2 > 1e-10 * lmax^2 || throw(ArgumentError(
        "shell element $(el.index) (:$(el.id)) is a degenerate triangle (collinear or coincident nodes)"))
    return nothing
end
