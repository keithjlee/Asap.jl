"""
    VariableElement{T} <: AbstractElement{T}

A **super-element**: one user-facing member whose cross-section varies
along its length as a chain of prismatic segments. The successor of the
legacy `BridgeElement`-era workflow of manually shattering members — here
the model is never mutated and no phantom nodes appear in `model.nodes`.

Interior joints between segments become **internal degrees of freedom**:
`process!` allocates them 6 slots each in the global DOF space *after* all
nodal slots. Assembly, solving, and recovery treat them exactly like nodal
DOFs (they are exact, not condensed — which keeps future mass and geometric
stiffness formulations exact too), but they are invisible to the user: the
element is queried as a single piece.

# Fields
- `nodeStart::Node{T}`, `nodeEnd::Node{T}`: the member's real end nodes;
  interior joints lie on the straight line between them
- `sections::Vector{<:AbstractSection{T}}`: one section per segment,
  ordered start → end
- `breaks::Vector{T}`: interior joint positions as strictly increasing
  fractions ∈ (0, 1) of the member length; `length(breaks) == length(sections) - 1`
- `ends::EndConditions{T}`: end conditions at the OUTER ends only (interior
  joints are rigid by construction)
- `Ψ::T`: roll angle, shared by all segments [rad]
- `id::Symbol`, `index::Int`: as for other elements
- `internal_offset::Int`: first global DOF slot of the interior-joint block
  (assigned by `process!`; internal bookkeeping)

# Constructors
    VariableElement(nodeStart, nodeEnd, sections, breaks, id = :variable;
                    release = :fixedfixed, Ψ = π/2)
    VariableElement(nodeStart, nodeEnd, sections, id; ...)   # equal segments

# Element-DOF layout
Slots 1–6: start node; 7–12: end node; 13 onward: interior joints in order.
[`segment_slots`](@ref) maps each segment's 12 DOFs into this layout.

# Examples
```julia-repl
julia> haunched = VariableElement(n1, n2, [deep, mid, shallow], [0.2, 0.5], :girder)

julia> stepped = VariableElement(n1, n2, [big, small], :column)  # break at 0.5
```
"""
mutable struct VariableElement{T} <: AbstractElement{T}
    nodeStart::Node{T}
    nodeEnd::Node{T}
    sections::Vector{AbstractSection{T}}
    breaks::Vector{T}
    ends::EndConditions{T}
    Ψ::T
    id::Symbol
    index::Int
    internal_offset::Int

    function VariableElement(nodeStart::Node{T}, nodeEnd::Node{T},
        sections::Vector{<:AbstractSection{T}}, breaks::Vector{<:Real},
        id::Symbol=:variable; release::Symbol=:fixedfixed, Ψ::Real=pi / 2) where {T}
        @assert length(sections) >= 2 "a VariableElement needs ≥ 2 segments (use FrameElement for prismatic members)"
        @assert length(breaks) == length(sections) - 1 "need one interior break per segment boundary"
        @assert issorted(breaks) && first(breaks) > 0 && last(breaks) < 1 "breaks must be strictly increasing in (0, 1)"
        return new{T}(nodeStart, nodeEnd, collect(AbstractSection{T}, sections),
            Vector{T}(breaks), EndConditions(release; T=T), T(Ψ), id, 0, 0)
    end
end

VariableElement(nodeStart::Node{T}, nodeEnd::Node{T},
    sections::Vector{<:AbstractSection{T}}, id::Symbol=:variable;
    kwargs...) where {T} =
    VariableElement(nodeStart, nodeEnd, sections,
        collect(range(0, 1; length=length(sections) + 1))[2:end-1], id; kwargs...)

nodes(el::VariableElement) = (el.nodeStart, el.nodeEnd)

"""
    n_segments(el) -> Int

Number of prismatic segments (1 for primitive elements).
"""
n_segments(::AbstractElement) = 1
n_segments(el::VariableElement) = length(el.sections)

n_internal_dofs(el::VariableElement) = 6 * (n_segments(el) - 1)

"""
    segment_fractions(el) -> Vector{T}

Segment boundary positions as fractions of the member length:
`[0, breaks..., 1]`.
"""
segment_fractions(el::VariableElement{T}) where {T} = vcat(zero(T), el.breaks, one(T))

"""
    segment_slots(el, s) -> Vector{Int}

The 12 element-DOF layout slots of segment `s` (its start block then end
block), in the element's slot layout (1–6 start node, 7–12 end node, 13+
interior joints).
"""
function segment_slots(el::VariableElement, s::Int)
    nseg = n_segments(el)
    startblock = s == 1 ? (1:6) : (12+6*(s-2)+1):(12+6*(s-1))
    endblock = s == nseg ? (7:12) : (12+6*(s-1)+1):(12+6*s)
    return vcat(collect(startblock), collect(endblock))
end

"""
    segment_ends(el, s) -> EndConditions

End conditions of segment `s`: the member's outer conditions at the outer
ends, rigid at interior joints.
"""
function segment_ends(el::VariableElement{T}, s::Int) where {T}
    r = rigid_end(T)
    return EndConditions(s == 1 ? el.ends.e1 : r,
        s == n_segments(el) ? el.ends.e2 : r)
end

"""
    dof_signature(el::VariableElement) -> NTuple

Blockwise activity over the element's full slot layout: translations always
active; outer rotation blocks follow the same blockwise logic as
`FrameElement` (using the outer end springs; the torsion chain releases if
either outer torsional spring is zero); interior joint blocks are fully
active (rigidly connected on both sides).
"""
function dof_signature(el::VariableElement)
    e1, e2 = el.ends.e1, el.ends.e2
    torsion = !(iszero(e1.kt) || iszero(e2.kt))
    rot1 = torsion || !iszero(e1.ky) || !iszero(e1.kz)
    rot2 = torsion || !iszero(e2.ky) || !iszero(e2.kz)
    outer = (true, true, true, rot1, rot1, rot1,
        true, true, true, rot2, rot2, rot2)
    return (outer..., ntuple(_ -> true, n_internal_dofs(el))...)
end

"""
    stiffness(el::VariableElement, sections, x1, x2) -> Matrix

Global stiffness over the element's full DOF layout: the sum of each
segment's [`frame_stiffness`](@ref) embedded at its [`segment_slots`](@ref).
Interior joint positions derive from the end positions (`xᵢ = x1 + tᵢ(x2−x1)`),
so geometry gradients flow through interior joints automatically. Built as
a pure sum of selector embeddings (AD-transparent); returns a dense Matrix
since the size varies with segment count.
"""
function stiffness(el::VariableElement, sections::AbstractVector,
    x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real})
    n = ndofs(el)
    ts = segment_fractions(el)
    v1 = SVector{3}(x1)
    v2 = SVector{3}(x2)

    parts = map(1:n_segments(el)) do s
        xa = v1 + ts[s] * (v2 - v1)
        xb = v1 + ts[s+1] * (v2 - v1)
        Ks = frame_stiffness(sections[s], segment_ends(el, s), xa, xb, el.Ψ)
        P = _segment_selector(n, segment_slots(el, s))
        P * Ks * P'
    end
    return sum(parts)
end

stiffness(el::VariableElement, x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real}) =
    stiffness(el, el.sections, x1, x2)

# n×12 selector embedding a segment's 12 DOFs into the element layout
function _segment_selector(n::Int, slots::Vector{Int})
    P = zeros(n, 12)
    for (j, s) in enumerate(slots)
        P[s, j] = 1.0
    end
    return P
end

"""
    segment_positions(el, s) -> (SVector{3}, SVector{3})

Current global positions of segment `s`'s two ends.
"""
function segment_positions(el::VariableElement, s::Int)
    ts = segment_fractions(el)
    v1 = el.nodeStart.position
    v2 = el.nodeEnd.position
    return (v1 + ts[s] * (v2 - v1), v1 + ts[s+1] * (v2 - v1))
end

"""
    locate_segment(el, t) -> (s, τ)

Map a fraction `t ∈ [0, 1]` of the whole member to its segment index `s`
and the local fraction `τ ∈ [0, 1]` within that segment — the resolution
step behind unified queries like `moment_at(el, 0.5)`.
"""
function locate_segment(el::VariableElement, t::Real)
    @assert 0 <= t <= 1 "position must be a fraction of the member length"
    ts = segment_fractions(el)
    s = min(searchsortedlast(ts, t) , n_segments(el))
    s = max(s, 1)
    τ = (t - ts[s]) / (ts[s+1] - ts[s])
    return s, τ
end

local_frame(el::VariableElement; tol::Real=1e-6) =
    local_frame(el.nodeStart.position, el.nodeEnd.position, el.Ψ; tol=tol)
