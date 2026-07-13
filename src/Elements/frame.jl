"""
    FrameElement{T, S<:AbstractSection{T}} <: AbstractElement{T}

A 3D frame (beam-column) element: carries axial force, biaxial bending,
transverse shear, and torsion between two nodes. The successor of the
legacy `Element` — renamed for what it is, and stripped of all analysis
state (stiffness, transformation, DOF ids, and forces live in the analysis
cache and results objects, not on the element).

Formulated as an Euler-Bernoulli member with **end springs**
([`EndConditions`](@ref)): classical releases are the spring limits, and
finite values model semi-rigid connections. Its stiffness couples all 12
nodal DOFs *except* those decoupled by zero end springs — reflected in its
[`dof_signature`](@ref), so fully released rotations never poison the
global system.

# Fields
- `nodeStart::Node{T}`, `nodeEnd::Node{T}`: connected nodes; the local x
  axis runs from start to end
- `section::S`: cross-section (any [`AbstractSection`](@ref)); mutable so
  section swaps in design iteration don't rebuild the element
- `ends::EndConditions{T}`: connection stiffnesses at both ends, in local
  coordinates
- `Ψ::T`: roll angle [rad] — rotation of the section about the element axis.
  Default `π/2` (legacy convention: local strong axis resists vertical load
  for typical horizontal members)
- `id::Symbol`: user tag for group queries
- `index::Int`: position in the model's element vector; set by `process!`

# Constructors
    FrameElement(nodeStart, nodeEnd, section, id = :element;
                 release = :fixedfixed, Ψ = π/2)
    FrameElement(nodeStart, nodeEnd, section, ends::EndConditions, id; Ψ = π/2)

# Examples
```julia-repl
julia> beam = FrameElement(n1, n2, sec, :girder)

julia> hinged = FrameElement(n1, n2, sec; release = :fixedfree)

julia> semirigid = FrameElement(n1, n2, sec,
           EndConditions(EndSprings(Inf, Inf, 5e4, 5e4), rigid_end()), :connection)
```
"""
mutable struct FrameElement{T,S<:AbstractSection{T}} <: AbstractElement{T}
    nodeStart::Node{T}
    nodeEnd::Node{T}
    section::S
    ends::EndConditions{T}
    Ψ::T
    id::Symbol
    index::Int

    function FrameElement(nodeStart::Node{T}, nodeEnd::Node{T}, section::S,
        ends::EndConditions, id::Symbol=:element; Ψ::Real=pi / 2) where {T,S<:AbstractSection{T}}
        e = EndConditions(EndSprings(T(ends.e1.kx), T(ends.e1.kt), T(ends.e1.ky), T(ends.e1.kz)),
            EndSprings(T(ends.e2.kx), T(ends.e2.kt), T(ends.e2.ky), T(ends.e2.kz)))
        return new{T,S}(nodeStart, nodeEnd, section, e, T(Ψ), id, 0)
    end
end

FrameElement(nodeStart::Node{T}, nodeEnd::Node{T}, section::AbstractSection{T},
    id::Symbol=:element; release::Symbol=:fixedfixed, Ψ::Real=pi / 2) where {T} =
    FrameElement(nodeStart, nodeEnd, section, EndConditions(release; T=T), id; Ψ=Ψ)

nodes(el::FrameElement) = (el.nodeStart, el.nodeEnd)

"""
    dof_signature(el::FrameElement) -> NTuple{12,Bool}

Which of the element's 12 nodal DOF slots its GLOBAL stiffness can couple.

Activity is decided **per rotation block, not per slot**: the local-to-global
rotation mixes torsion and bending rotations within a node's 3-slot rotation
block, so a single released local rotation still couples all three global
rotational DOFs through the other two. A node's rotation block is inactive
only when the element's entire local rotational row block at that end is
zero — i.e. both bending springs released AND torsion released (torsion
releases whenever either end's torsional spring is zero, since a torsion
chain with a free end restrains nothing).

Consequences: a `:fixedfree` element leaves its far node's rotations
entirely untouched (no singular mode if nothing else connects there —
Keith's hinge/torsion pain point, solved structurally); a `:joist` element
keeps torsional coupling and therefore marks whole rotation blocks active.
Translations are always coupled.
"""
function dof_signature(el::FrameElement)
    e1, e2 = el.ends.e1, el.ends.e2
    torsion = !(iszero(e1.kt) || iszero(e2.kt))
    rot1 = torsion || !iszero(e1.ky) || !iszero(e1.kz)
    rot2 = torsion || !iszero(e2.ky) || !iszero(e2.kz)
    return (
        true, true, true, rot1, rot1, rot1,
        true, true, true, rot2, rot2, rot2,
    )
end

"""
    stiffness(el::FrameElement, x1, x2) -> SMatrix{12,12}

Global-coordinate stiffness of the frame element with its end positions
passed explicitly (the fast path passes node positions; the AD path passes
entries of a differentiable state). Delegates to the pure
[`frame_stiffness`](@ref) kernel.
"""
stiffness(el::FrameElement, x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real}) =
    frame_stiffness(el.section, el.ends, x1, x2, el.Ψ)

"""
    TrussElement{T, S<:AbstractSection{T}} <: AbstractElement{T}

An axial-only (two-force) element: resists elongation along its axis and
nothing else. Its [`dof_signature`](@ref) touches only the six nodal
translations — rotational DOF slots are simply never marked active, which
is how truss-only models solve a translations-only system and how truss and
frame elements mix freely in one model.

# Fields
- `nodeStart::Node{T}`, `nodeEnd::Node{T}`: connected nodes
- `section::S`: cross-section — only `EA` and `ρA` are ever queried, so an
  axial-only `Section(material, A)` suffices
- `id::Symbol`: user tag
- `index::Int`: set by `process!`

# Constructor
    TrussElement(nodeStart, nodeEnd, section, id = :element)
"""
mutable struct TrussElement{T,S<:AbstractSection{T}} <: AbstractElement{T}
    nodeStart::Node{T}
    nodeEnd::Node{T}
    section::S
    id::Symbol
    index::Int

    function TrussElement(nodeStart::Node{T}, nodeEnd::Node{T}, section::S,
        id::Symbol=:element) where {T,S<:AbstractSection{T}}
        return new{T,S}(nodeStart, nodeEnd, section, id, 0)
    end
end

nodes(el::TrussElement) = (el.nodeStart, el.nodeEnd)

"""
    dof_signature(el::TrussElement) -> NTuple{12,Bool}

Truss elements couple only the translational slots of their two nodes;
all six rotational slots are `false` and never activate global DOFs.
"""
dof_signature(::TrussElement) = (true, true, true, false, false, false,
    true, true, true, false, false, false)

"""
    stiffness(el::TrussElement, x1, x2) -> SMatrix{12,12}

Global-coordinate stiffness of the truss element, embedded in the full
12-slot (two-node) DOF layout with zeros on all rotational slots. The
active 6×6 translational block is the pure [`truss_stiffness`](@ref)
kernel; assembly consults [`dof_signature`](@ref) so the zero rotational
rows are never scattered into the global matrix.
"""
function stiffness(el::TrussElement, x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real})
    B = truss_stiffness(el.section, x1, x2)          # 6×6 on translations
    T = eltype(B)
    k = _mutable12(T)
    tidx = SVector(1, 2, 3, 7, 8, 9)                  # translational slots
    @inbounds for a in 1:6, b in 1:6
        k[tidx[a], tidx[b]] = B[a, b]
    end
    return SMatrix{12,12,T}(k)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shared element utilities
# ─────────────────────────────────────────────────────────────────────────────

"""
    endpoints(el) -> (SVector{3}, SVector{3})

The element's start and end positions.
"""
endpoints(el::AbstractElement) = (el.nodeStart.position, el.nodeEnd.position)

"""
    midpoint(el) -> SVector{3}

Position of the element's midpoint.
"""
midpoint(el::AbstractElement) = (el.nodeStart.position + el.nodeEnd.position) / 2

"""
    length(el) -> T

Length of the element [length].
"""
Base.length(el::AbstractElement) = element_length(el.nodeStart.position, el.nodeEnd.position)

"""
    local_frame(el; tol = 1e-6) -> SMatrix{3,3}

The element's local coordinate frame (rows = local x, y, z in global
coordinates). Frame elements use their roll angle `Ψ`; truss elements have
no roll (Ψ = 0).
"""
local_frame(el::FrameElement; tol::Real=1e-6) =
    local_frame(el.nodeStart.position, el.nodeEnd.position, el.Ψ; tol=tol)
local_frame(el::TrussElement; tol::Real=1e-6) =
    local_frame(el.nodeStart.position, el.nodeEnd.position, 0.0; tol=tol)
