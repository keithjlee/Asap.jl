"""
    AbstractLoad{T<:Real}

Supertype of all loads. Loads are immutable, reference their target (node
or element) directly, and carry a `case::Symbol` tag grouping them into a
load case for later combination (default `:LC1`).

Two families:
- [`NodeLoad`](@ref): concentrated actions applied directly to nodal DOFs
- [`ElementLoad`](@ref): actions applied along an element, converted to
  equivalent (fixed-end) nodal forces during load assembly

The single extension point for new element load types is
[`fixed_end_forces`](@ref)`(load, section, ends, x1, x2, rollangle)` returning the
clamped-clamped local 12-vector — end-condition condensation and rotation to
global coordinates are applied generically afterwards.
"""
abstract type AbstractLoad{T<:Real} end

"Concentrated actions applied directly to a node's DOFs. See [`AbstractLoad`](@ref)."
abstract type NodeLoad{T} <: AbstractLoad{T} end

"Actions applied along an element, lowered to fixed-end forces. See [`AbstractLoad`](@ref)."
abstract type ElementLoad{T} <: AbstractLoad{T} end

"""
    NodeForce{T} <: NodeLoad{T}

A concentrated force `[Fx, Fy, Fz]` [force] applied to a node in GLOBAL
coordinates.

    NodeForce(node, value; id = :force, case = :LC1)
"""
struct NodeForce{T} <: NodeLoad{T}
    node::Node{T}
    value::SVector{3,T}
    id::Symbol
    case::Symbol

    function NodeForce(node::Node{T}, value::AbstractVector{<:Real};
        id::Symbol=:force, case::Symbol=:LC1) where {T}
        @assert length(value) == 3 "force vector must be [Fx, Fy, Fz] in global coordinates"
        return new{T}(node, SVector{3,T}(value), id, case)
    end
end

"""
    NodeMoment{T} <: NodeLoad{T}

A concentrated moment `[Mx, My, Mz]` [force·length] applied to a node's
rotational DOFs in GLOBAL coordinates.

    NodeMoment(node, value; id = :moment, case = :LC1)
"""
struct NodeMoment{T} <: NodeLoad{T}
    node::Node{T}
    value::SVector{3,T}
    id::Symbol
    case::Symbol

    function NodeMoment(node::Node{T}, value::AbstractVector{<:Real};
        id::Symbol=:moment, case::Symbol=:LC1) where {T}
        @assert length(value) == 3 "moment vector must be [Mx, My, Mz] in global coordinates"
        return new{T}(node, SVector{3,T}(value), id, case)
    end
end

"""
    DistributedLoad{T} <: ElementLoad{T}

The canonical distributed load: a **piecewise-linear intensity** along the
element, in a fixed direction. Every distributed loading — uniform,
trapezoidal, triangular, partial-span, tributary — lowers to this one type,
and one fixed-end-force integration handles them all.

# Fields
- `element`: the loaded element
- `t::Vector{T}`: sorted breakpoint positions as fractions ∈ [0, 1] of the
  element length. Intensity is zero outside `[t[1], t[end]]`, so partial-span
  loads need no special casing.
- `w::Vector{T}`: load intensity at each breakpoint [force/length], varying
  linearly between breakpoints
- `direction::SVector{3,T}`: unit direction of the load
- `coords::Symbol`: `:global` (direction fixed in space — gravity, wind) or
  `:local` (direction follows the element's local axes)
- `id`, `case`: tags

# Convenience constructors
    LineLoad(element, wvec; id, case)              # uniform, full span, global
    DistributedLoad(element, t, w, direction; coords = :global, id, case)

# Examples
```julia-repl
julia> LineLoad(beam, [0.0, 0.0, -2.0])                       # uniform gravity-direction

julia> DistributedLoad(beam, [0.2, 0.8], [0.0, 5.0], SVector(0.0, 0.0, -1.0))
       # partial-span triangular ramp between 20% and 80% of the span
```
"""
struct DistributedLoad{T} <: ElementLoad{T}
    element::AbstractElement{T}
    t::Vector{T}
    w::Vector{T}
    direction::SVector{3,T}
    coords::Symbol
    id::Symbol
    case::Symbol

    function DistributedLoad(element::AbstractElement{T}, t::AbstractVector{<:Real},
        w::AbstractVector{<:Real}, direction::AbstractVector{<:Real};
        coords::Symbol=:global, id::Symbol=:distributed, case::Symbol=:LC1) where {T}
        @assert length(t) == length(w) >= 2 "need matching breakpoints t and intensities w (≥ 2)"
        @assert issorted(t) "breakpoints t must be sorted"
        @assert first(t) >= 0 && last(t) <= 1 "breakpoints must lie in [0, 1]"
        @assert coords in (:global, :local) "coords must be :global or :local"
        d = SVector{3,T}(normalize(direction))
        return new{T}(element, Vector{T}(t), Vector{T}(w), d, coords, id, case)
    end
end

"""
    LineLoad(element, wvec; id = :lineload, case = :LC1)

A uniform full-span distributed load given as a GLOBAL vector `[wx, wy, wz]`
[force/length] — the legacy `LineLoad` signature, lowered to a
[`DistributedLoad`](@ref) with magnitude `|wvec|` along `normalize(wvec)`.
"""
function LineLoad(element::AbstractElement{T}, wvec::AbstractVector{<:Real};
    id::Symbol=:lineload, case::Symbol=:LC1) where {T}
    @assert length(wvec) == 3 "load vector must be [wx, wy, wz] in global coordinates"
    mag = norm(wvec)
    @assert mag > 0 "zero line load"
    return DistributedLoad(element, T[0, 1], [T(mag), T(mag)], wvec ./ mag;
        coords=:global, id=id, case=case)
end

"""
    TrapezoidLoad(element, t1, t2, w1, w2, direction;
                  coords = :global, id = :trapezoidload, case = :LC1)

A partial-span trapezoidal distributed load: intensity varying linearly
from `w1` [force/length] at fraction `t1` to `w2` at fraction `t2` of the
element length, zero elsewhere, acting along the unit `direction`. A
triangular load is the `w1 = 0` (or `w2 = 0`) special case.

Lowered to the canonical [`DistributedLoad`](@ref) — one integration engine
covers every distributed shape.

# Examples
```julia-repl
julia> TrapezoidLoad(beam, 0.0, 1.0, 2.0, 5.0, [0.0, 0.0, -1.0])   # full-span trapezoid

julia> TrapezoidLoad(beam, 0.3, 0.9, 0.0, 4.0, [0.0, -1.0, 0.0])   # partial triangular ramp
```
"""
function TrapezoidLoad(element::AbstractElement{T}, t1::Real, t2::Real,
    w1::Real, w2::Real, direction::AbstractVector{<:Real};
    coords::Symbol=:global, id::Symbol=:trapezoidload, case::Symbol=:LC1) where {T}
    @assert t1 < t2 "need t1 < t2 for a trapezoid span"
    return DistributedLoad(element, T[t1, t2], T[w1, w2], direction;
        coords=coords, id=id, case=case)
end

"""
    PointMoment{T} <: ElementLoad{T}

A concentrated moment `[Mx, My, Mz]` [force·length] applied along an
element at fraction `position` ∈ (0, 1) of its length — a capability the
legacy library lacked. Its fixed-end forces pair the moment with the
SLOPES of the shape functions (a moment does work through rotation), which
is why the FEF kernel needs the Hermite derivatives.

    PointMoment(element, position, value; coords = :global, id = :pointmoment, case = :LC1)
"""
struct PointMoment{T} <: ElementLoad{T}
    element::AbstractElement{T}
    position::T
    value::SVector{3,T}
    coords::Symbol
    id::Symbol
    case::Symbol

    function PointMoment(element::AbstractElement{T}, position::Real,
        value::AbstractVector{<:Real};
        coords::Symbol=:global, id::Symbol=:pointmoment, case::Symbol=:LC1) where {T}
        @assert 0 < position < 1 "position must be ∈ (0, 1) — a fraction of the element length"
        @assert length(value) == 3 "moment vector must be [Mx, My, Mz]"
        @assert coords in (:global, :local) "coords must be :global or :local"
        return new{T}(element, T(position), SVector{3,T}(value), coords, id, case)
    end
end

"""
    PointLoad{T} <: ElementLoad{T}

A concentrated force `[Px, Py, Pz]` applied along an element at fraction
`position` ∈ (0, 1) of its length.

    PointLoad(element, position, value; coords = :global, id = :pointload, case = :LC1)
"""
struct PointLoad{T} <: ElementLoad{T}
    element::AbstractElement{T}
    position::T
    value::SVector{3,T}
    coords::Symbol
    id::Symbol
    case::Symbol

    function PointLoad(element::AbstractElement{T}, position::Real,
        value::AbstractVector{<:Real};
        coords::Symbol=:global, id::Symbol=:pointload, case::Symbol=:LC1) where {T}
        @assert 0 < position < 1 "position must be ∈ (0, 1) — a fraction of the element length"
        @assert length(value) == 3 "load vector must be [Px, Py, Pz]"
        @assert coords in (:global, :local) "coords must be :global or :local"
        return new{T}(element, T(position), SVector{3,T}(value), coords, id, case)
    end
end

"""
    SelfWeight{T} <: ElementLoad{T}

Self-weight of an element: a distributed load of intensity
`ρA(section) · |g|` in the direction of `g`, where `ρA` is the section's
mass per unit length. There is exactly one source of truth for mass — the
section accessor — so density lives on the `Material`/section, never on the
load.

    SelfWeight(element; g = SVector(0, 0, -9.80665), factor = 1, id, case)

`g` is the gravitational acceleration vector [length/time²] in global
coordinates (default: SI standard gravity in −Z); `factor` is a plain
multiplier (leave load-combination factors to `LoadCombination`s).
"""
struct SelfWeight{T} <: ElementLoad{T}
    element::AbstractElement{T}
    g::SVector{3,T}
    factor::T
    id::Symbol
    case::Symbol

    function SelfWeight(element::AbstractElement{T};
        g::AbstractVector{<:Real}=SVector(0.0, 0.0, -9.80665),
        factor::Real=1, id::Symbol=:selfweight, case::Symbol=:LC1) where {T}
        @assert length(g) == 3 "g must be a 3-vector of gravitational acceleration"
        return new{T}(element, SVector{3,T}(g), T(factor), id, case)
    end
end
