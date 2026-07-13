"""
    EndSprings{T<:Real}

Connection stiffness at ONE end of a frame element, expressed in the
element's local coordinate system.

This generalizes the classical binary release: instead of a DOF being either
rigidly tied to the node or fully released, each end carries a finite spring
stiffness for the four end actions that can meaningfully be released or
softened in a line element:

- `kx`: axial spring [force/length] — stiffness of the axial connection
- `kt`: torsional spring [force·length/rad] — twist about the element axis
- `ky`: rotational spring about the local y axis [force·length/rad] —
  bending in the local x–z plane
- `kz`: rotational spring about the local z axis [force·length/rad] —
  bending in the local x–y plane

The limits recover classical behavior exactly: `Inf` = rigid connection
(the DOF is fully coupled to the node), `0` = ideal release (hinge/slide).
Intermediate values model **semi-rigid connections** — bolted end plates,
concrete joint regions with partial fixity, etc. — and are legitimate,
differentiable design variables.

Transverse shear releases are deliberately not supported (as in the legacy
library): elastic supports are modeled with nodal springs instead.

See also [`EndConditions`](@ref), [`RELEASES`](@ref).
"""
struct EndSprings{T<:Real}
    kx::T # axial
    kt::T # torsion
    ky::T # rotation about local y
    kz::T # rotation about local z

    function EndSprings(kx::Real, kt::Real, ky::Real, kz::Real)
        args = promote(kx, kt, ky, kz)
        @assert all(k -> k >= 0, args) "spring stiffnesses must be ≥ 0 (use Inf for rigid)"
        return new{eltype(args)}(args...)
    end
end

"""
    rigid_end(T = Float64) -> EndSprings{T}

The rigid connection: all four end springs infinite — the classical fully
fixed element end.
"""
rigid_end(::Type{T}=Float64) where {T<:Real} = EndSprings(T(Inf), T(Inf), T(Inf), T(Inf))

"""
    pinned_end(T = Float64) -> EndSprings{T}

The classical hinge: axial rigid, all rotational stiffnesses zero (torsion
and both bending rotations released).
"""
pinned_end(::Type{T}=Float64) where {T<:Real} = EndSprings(T(Inf), zero(T), zero(T), zero(T))

"""
    EndConditions{T<:Real}

The pair of [`EndSprings`](@ref) at the start (`e1`) and end (`e2`) of a
frame element, in element local coordinates.

Constructed either explicitly from two `EndSprings`, or from a classical
release symbol (see [`RELEASES`](@ref)):

    EndConditions(e1::EndSprings, e2::EndSprings)
    EndConditions(release::Symbol; T = Float64)

# Examples
```julia-repl
julia> EndConditions(:fixedfree)          # hinge at the far end

julia> EndConditions(EndSprings(Inf, Inf, 5e4, 5e4),   # semi-rigid start
                     rigid_end())                      # rigid end
```
"""
struct EndConditions{T<:Real}
    e1::EndSprings{T} # start-node end
    e2::EndSprings{T} # end-node end
end

function EndConditions(e1::EndSprings{T1}, e2::EndSprings{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    EndConditions{T}(EndSprings(T(e1.kx), T(e1.kt), T(e1.ky), T(e1.kz)),
        EndSprings(T(e2.kx), T(e2.kt), T(e2.ky), T(e2.kz)))
end

function EndConditions(release::Symbol; T::Type{<:Real}=Float64)
    @assert haskey(RELEASES, release) "unknown release :$release — choose from $(sort!(collect(keys(RELEASES))))"
    e1, e2 = RELEASES[release]
    make = e -> EndSprings(T(e[1]), T(e[2]), T(e[3]), T(e[4]))
    return EndConditions(make(e1), make(e2))
end

"""
    RELEASES

Map from the classical release symbols to their exact end-spring limits,
as `(start, end)` tuples of `(kx, kt, ky, kz)`. These are the same five
release types the legacy library implemented as distinct stiffness
matrices; here they are just data:

| symbol | start end | far end | classical meaning |
|---|---|---|---|
| `:fixedfixed` | rigid | rigid | fully continuous member (default) |
| `:fixedfree`  | rigid | rotations released | hinge at far end |
| `:freefixed`  | rotations released | rigid | hinge at start end |
| `:freefree`   | rotations released | rotations released | axial-only (truss-like) |
| `:joist`      | bending released, torsion kept | same | joist/purlin idealization |
"""
const RELEASES = Dict{Symbol,NTuple{2,NTuple{4,Float64}}}(
    :fixedfixed => ((Inf, Inf, Inf, Inf), (Inf, Inf, Inf, Inf)),
    :fixedfree => ((Inf, Inf, Inf, Inf), (Inf, 0.0, 0.0, 0.0)),
    :freefixed => ((Inf, 0.0, 0.0, 0.0), (Inf, Inf, Inf, Inf)),
    :freefree => ((Inf, 0.0, 0.0, 0.0), (Inf, 0.0, 0.0, 0.0)),
    :joist => ((Inf, Inf, 0.0, 0.0), (Inf, Inf, 0.0, 0.0)),
)

"""
    release_symbol(ends::EndConditions) -> Union{Symbol, Nothing}

If the end conditions exactly match one of the classical [`RELEASES`](@ref),
return its symbol; otherwise return `nothing` (a genuinely semi-rigid
connection).
"""
function release_symbol(ends::EndConditions)
    tup = e -> (e.kx, e.kt, e.ky, e.kz)
    for (sym, (r1, r2)) in RELEASES
        if tup(ends.e1) == r1 && tup(ends.e2) == r2
            return sym
        end
    end
    return nothing
end
