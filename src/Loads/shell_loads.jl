"""
    AreaLoad{T} <: ElementLoad{T}

A uniform surface traction on a [`ShellElement`](@ref) [force/length²] —
pressure on a slab or wall, wind on a façade panel, superimposed dead load.

Lowered to CONSISTENT nodal loads from the element's own interpolation
(including the work done on the drilling rotations by in-plane traction, and
the condensed internal modes); see `shell_load_vector`.

# Fields
- `element::AbstractElement{T}`: the loaded shell element
- `q::SVector{3,T}`: traction [force/length²], in global or local axes
- `frame::Symbol`: `:global` (default) or `:local` (element axes `e₁, e₂, e₃`;
  a local `q = (0, 0, p)` is a pressure `p` along the element normal and
  follows the element's orientation)
- `id::Symbol`, `case::Symbol`: tags

# Constructors
    AreaLoad(element, q::AbstractVector; frame = :global, id = :area, case = :LC1)
    AreaLoad(element, p::Real; id = :area, case = :LC1)   # normal pressure along +e₃

# Examples
```julia-repl
julia> AreaLoad(slab, [0.0, 0.0, -0.1])        # 0.1 [force/length²] downward
julia> AreaLoad(facade, -0.03; case = :W)       # suction along −e₃
```

`SelfWeight(element)` also applies to shells: `q = ρt·g·factor` from the
section's [`areal_mass`](@ref).
"""
struct AreaLoad{T} <: ElementLoad{T}
    element::AbstractElement{T}
    q::SVector{3,T}
    frame::Symbol
    id::Symbol
    case::Symbol

    function AreaLoad(element::AbstractElement{T}, q::AbstractVector{<:Real};
        frame::Symbol = :global, id::Symbol = :area, case::Symbol = :LC1) where {T}
        element isa ShellElement || throw(ArgumentError(
            "AreaLoad applies to ShellElements (got a $(nameof(typeof(element))))"))
        length(q) == 3 || throw(ArgumentError("AreaLoad traction must be a 3-vector"))
        frame in (:global, :local) || throw(ArgumentError("frame must be :global or :local"))
        return new{T}(element, SVector{3,T}(q), frame, id, case)
    end
end

AreaLoad(element::AbstractElement{T}, p::Real; id::Symbol = :area, case::Symbol = :LC1) where {T} =
    AreaLoad(element, SVector{3,T}(0, 0, p); frame = :local, id = id, case = case)

# the load's uniform traction in GLOBAL axes on its shell element
_shell_traction(load::AreaLoad, el::ShellElement) =
    load.frame === :global ? load.q : local_frame(el)' * load.q
_shell_traction(load::SelfWeight, el::ShellElement) =
    areal_mass(el.section) * load.factor * load.g
_shell_traction(load::ElementLoad, el::ShellElement) = throw(ArgumentError(
    "$(nameof(typeof(load))) cannot be applied to a ShellElement — use AreaLoad, SelfWeight or NodeForce"))
