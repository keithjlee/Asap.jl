"""
    LinkElement{T} <: AbstractElement{T}

A zero-length 2-node spring connecting two COINCIDENT nodes with six
uncoupled stiffnesses in a chosen frame — the connection element for precast
panel joints (horizontal/vertical joints of shear walls, bearing pads, dowel
connections), base isolators, or any "these two parts are attached, but not
rigidly" condition.

Element stiffness (global): `[k −k; −k k]` with
`k = blockdiag(Rᵀ·diag(kₜ)·R, Rᵀ·diag(kᵣ)·R)`, where `R` (rows = link axes)
orients the springs. Zero length is enforced at `process!`: a spring between
OFFSET nodes would need rigid arms for moment equilibrium; model those
explicitly (a stiff `FrameElement`) if you need them.

# Fields
- `nodeStart::Node{T}`, `nodeEnd::Node{T}`: the two coincident nodes
- `k::SVector{6,T}`: stiffnesses `(k₁, k₂, k₃, kθ₁, kθ₂, kθ₃)` along/about the
  link axes [force/length] and [force·length/rad]; 0 = released
- `R::SMatrix{3,3,T}`: rows are the link axes in global coordinates
  (default: global X, Y, Z)
- `id::Symbol`, `index::Int`

# Constructors
    LinkElement(n1, n2, k::AbstractVector; axes = I, id = :link)

# Examples
```julia-repl
julia> joint = LinkElement(top_of_lower_panel, bottom_of_upper_panel,
                           [1e6, 1e6, 5e6, 0, 0, 0]; id = :hjoint)   # shear + bearing, no moment
```
"""
mutable struct LinkElement{T} <: AbstractElement{T}
    nodeStart::Node{T}
    nodeEnd::Node{T}
    k::SVector{6,T}
    R::SMatrix{3,3,T,9}
    id::Symbol
    index::Int

    function LinkElement(n1::Node{T}, n2::Node{T}, k::AbstractVector{<:Real};
        axes::AbstractMatrix = SMatrix{3,3,T}(I), id::Symbol = :link) where {T}
        length(k) == 6 || throw(ArgumentError("LinkElement needs 6 stiffnesses"))
        all(>=(0), k) || throw(ArgumentError("LinkElement stiffnesses must be non-negative"))
        any(>(0), k) || throw(ArgumentError("LinkElement with all-zero stiffness connects nothing"))
        R = SMatrix{3,3,T}(axes)
        R * R' ≈ I || throw(ArgumentError("LinkElement axes must be orthonormal (rows)"))
        n1 === n2 && throw(ArgumentError("LinkElement needs two distinct nodes"))
        return new{T}(n1, n2, SVector{6,T}(k), R, id, 0)
    end
end

nodes(el::LinkElement) = (el.nodeStart, el.nodeEnd)

"""
    dof_signature(el::LinkElement) -> NTuple{12, Bool}

Activity per 3-DOF block (translations, rotations) of each node: a block is
active when ANY of its link stiffnesses is nonzero (the link axes may mix
global components, so activity cannot be decided per component).
"""
function dof_signature(el::LinkElement)
    t = any(>(0), el.k[SOneTo(3)])
    r = any(>(0), el.k[SVector(4, 5, 6)])
    return (t, t, t, r, r, r, t, t, t, r, r, r)
end

"""
    stiffness(el::LinkElement, x1, x2) -> SMatrix{12,12}

`[k −k; −k k]` in global coordinates (independent of position: zero length).
"""
function stiffness(el::LinkElement{T}, x1, x2) where {T}
    R = el.R
    kt = R' * Diagonal(el.k[SOneTo(3)]) * R
    kr = R' * Diagonal(el.k[SVector(4, 5, 6)]) * R
    Z = zero(SMatrix{3,3,T})
    k6 = vcat(hcat(kt, Z), hcat(Z, kr))
    return vcat(hcat(k6, -k6), hcat(-k6, k6))
end

function validate_rigidities(el::LinkElement)
    d = norm(el.nodeEnd.position - el.nodeStart.position)
    scale = max(norm(el.nodeStart.position), one(d))
    d <= 1e-9 * scale || throw(ArgumentError(
        "link element $(el.index) (:$(el.id)) connects non-coincident nodes (distance $d); " *
        "LinkElement is zero-length — use a FrameElement for offset connections"))
    return nothing
end
