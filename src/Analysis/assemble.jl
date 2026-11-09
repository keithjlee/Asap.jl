#=
Numeric assembly — the in-place fast path.

The sparsity pattern was frozen by the symbolic pass; assembling a new
stiffness state is: zero nzval, evaluate each element's stiffness kernel
(stack SMatrix), scatter through the precomputed nzmap, add spring
diagonals. No COO, no sparse() rebuild, no allocation.
=#

"""
    assemble_K!(cache::AnalysisCache, model) -> cache.K

Numerically (re)assemble the free×free global stiffness matrix in place.
Element kernels read current node positions and section properties, so
geometry/section changes since the last call are picked up automatically —
only *topology* changes require re-processing.
"""
function assemble_K!(cache::AnalysisCache{T}, model::Model{T}) where {T}
    nz = nonzeros(cache.K)
    fill!(nz, zero(T))
    for g in cache.groups
        _assemble_group!(nz, g)               # function barrier: concrete E
    end
    @inbounds for (pos, si, i) in cache.spring_entries   # values read FRESH
        nz[pos] += T(model.springs[si].stiffness[i])
    end
    return cache.K
end

"""
    _refresh_spring_nzvec!(cache, model)

Rebuild the pure path's spring-diagonal vector from the model's CURRENT
springs (called by [`extract_state`](@ref) so replaced spring values are
picked up — mirrors the fresh read the in-place assembler does).
"""
function _refresh_spring_nzvec!(cache::AnalysisCache{T}, model::Model{T}) where {T}
    fill!(cache.spring_nzvec, zero(T))
    @inbounds for (pos, si, i) in cache.spring_entries
        cache.spring_nzvec[pos] += T(model.springs[si].stiffness[i])
    end
    return cache.spring_nzvec
end

function _assemble_group!(nz::Vector{T}, g::ElementGroup{E}) where {T,E}
    @inbounds for (e, el) in enumerate(g.elements)
        Ke = stiffness(el, el.nodeStart.position, el.nodeEnd.position)
        slots = g.slots[e]
        m = g.nzmap[e]
        n = length(slots)
        p = 0
        for b in 1:n
            sb = slots[b]
            for a in 1:n
                p += 1
                pos = m[p]
                pos > 0 && (nz[pos] += Ke[slots[a], sb])
            end
        end
    end
    return nz
end

"""
    assemble_loads!(cache::AnalysisCache, model) -> (cache.P, cache.Pf)

Build the full-space nodal load vector `P` and fixed-end force vector `Pf`
from the model's loads.

Node loads accumulate directly into `P` (erroring usefully if a component
targets an *inactive* DOF — e.g. a moment on a node connected only to truss
elements — instead of surfacing later as a bare `SingularException`).

Element loads run the generic lowering: clamped local FEFs from the load's
[`fixed_end_forces`](@ref) kernel → end-condition condensation
([`condense_fef`](@ref)) → accumulate locally in `cache.q_local` (for force
recovery) and globally (blockwise Λᵀ rotation) into `Pf`.

Also rebuilds `cache.fef` ([`FEFTable`](@ref)), the plain-data mirror of the
applied element loads that the pure path re-evaluates per design state.
"""
function assemble_loads!(cache::AnalysisCache{T}, model::Model{T};
    case::Union{Nothing,Symbol}=nothing) where {T}
    fill!(cache.P, zero(T))
    fill!(cache.Pf, zero(T))
    for qs in cache.q_local, q in qs
        fill!(q, zero(T))
    end

    inactive = Set(cache.partition.inactive)
    tracked = Any[]
    static = zeros(T, length(cache.Pf))

    for load in model.loads
        case === nothing || load.case === case || continue
        if load isa ElementLoad && load.element isa VariableElement
            before = copy(cache.Pf)
            _apply_load!(cache, load, inactive)
            static .+= cache.Pf .- before
        else
            _apply_load!(cache, load, inactive)
            load isa ElementLoad && push!(tracked, load)
        end
    end
    cache.fef = _build_fef_table(T, tracked, length(cache.Pf), static)
    return cache.P, cache.Pf
end

# the plain-data table of state-tracked element loads (FrameElement only)
function _build_fef_table(::Type{T}, loads::Vector{Any}, n_global::Int, static::Vector{T}) where {T}
    n = length(loads)
    elem = Vector{Int}(undef, n); i1 = similar(elem); i2 = similar(elem)
    roll = Vector{T}(undef, n); ends = Vector{Any}(undef, n)
    I_ = Int[]; J_ = Int[]
    for (j, load) in enumerate(loads)
        el = load.element::FrameElement
        elem[j] = el.index
        i1[j] = el.nodeStart.index
        i2[j] = el.nodeEnd.index
        roll[j] = el.rollangle
        ends[j] = el.ends
        g = element_global_dofs(el)
        for k in 1:12
            push!(I_, g[k]); push!(J_, 12(j - 1) + k)
        end
    end
    scatter = sparse(I_, J_, ones(T, length(I_)), n_global, 12n)
    return FEFTable{T}(loads, elem, i1, i2, roll, ends, scatter, static)
end

"""
    condensed_fef(load, section, ends, x1, x2, rollangle) -> SVector{12}

The element load's fixed-end force vector in LOCAL coordinates, consistent
with the element's end conditions: the clamped kernel
([`fixed_end_forces`](@ref)) followed by end-spring condensation
([`condense_fef`](@ref)). This is `q̃` in the local force recovery
`f = k·u_local + q̃`. Pure and differentiable in the positions `x1`, `x2`,
the `section` and the `ends`.
"""
function condensed_fef(load::ElementLoad, section::AbstractSection, ends::EndConditions,
    x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real}, rollangle::Real)
    q = fixed_end_forces(load, section, ends, x1, x2, rollangle)
    return condense_fef(q, section, element_length(x1, x2), ends)
end

"""
    element_fef_global(load, section, ends, x1, x2, rollangle) -> SVector{12}

[`condensed_fef`](@ref) rotated to GLOBAL coordinates (blockwise Λᵀ): the
element load's contribution to `Pf` at the element's 12 global DOFs. Pure
and differentiable; the stateful assembler and the pure path
([`assemble_Pf`](@ref)) both evaluate exactly this.
"""
function element_fef_global(load::ElementLoad, section::AbstractSection, ends::EndConditions,
    x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real}, rollangle::Real)
    q̃ = condensed_fef(load, section, ends, x1, x2, rollangle)
    return _rotate12(local_frame(x1, x2, rollangle)', q̃)
end

function _apply_load!(cache::AnalysisCache, load::NodeForce, inactive)
    s = node_dof_start(load.node)
    for i in 1:3
        g = s + i
        if !iszero(load.value[i]) && g in inactive
            error("NodeForce :$(load.id) applies a force to inactive DOF $(("Tx","Ty","Tz")[i]) " *
                  "of node :$(load.node.id) — no element or spring provides stiffness there")
        end
        cache.P[g] += load.value[i]
    end
end

function _apply_load!(cache::AnalysisCache, load::NodeMoment, inactive)
    s = node_dof_start(load.node)
    for i in 1:3
        g = s + 3 + i
        if !iszero(load.value[i]) && g in inactive
            error("NodeMoment :$(load.id) applies a moment about $(("X","Y","Z")[i]) " *
                  "to node :$(load.node.id), but that rotational DOF is inactive " *
                  "(connected elements provide no rotational stiffness — e.g. truss members or released ends)")
        end
        cache.P[g] += load.value[i]
    end
end

function _apply_load!(cache::AnalysisCache{T}, load::ElementLoad{T}, inactive) where {T}
    el = load.element
    el isa Union{FrameElement,VariableElement} ||
        error("element loads require a FrameElement or VariableElement " *
              "(a $(nameof(typeof(el))) cannot carry transverse load) — " *
              "apply NodeForces or subdivide instead")
    _apply_element_load!(cache, load, el)
end

function _apply_element_load!(cache::AnalysisCache{T}, load::ElementLoad{T},
    el::FrameElement) where {T}
    x1 = el.nodeStart.position
    x2 = el.nodeEnd.position

    # the same pure kernels the differentiable path evaluates (assemble_Pf)
    q̃ = condensed_fef(load, el.section, el.ends, x1, x2, el.rollangle)
    cache.q_local[el.index][1] .+= q̃

    Qg = _rotate12(local_frame(x1, x2, el.rollangle)', q̃)   # local → global, blockwise
    g = element_global_dofs(el)
    @inbounds for i in 1:12
        cache.Pf[g[i]] += Qg[i]
    end
end

"""
Element loads on a VariableElement are split across its segments — each
segment gets the load restricted to its span, in segment-local fractions —
then lowered per segment exactly like a primitive element, scattering into
the outer-node and internal-joint global DOFs.
"""
function _apply_element_load!(cache::AnalysisCache{T}, load::ElementLoad{T},
    el::VariableElement) where {T}
    g = element_global_dofs(el)
    ts = segment_fractions(el)
    Λ = local_frame(el)

    for s in 1:n_segments(el)
        subload = _restrict_load(load, ts[s], ts[s+1], el)
        subload === nothing && continue

        xa, xb = segment_positions(el, s)
        Ls = element_length(xa, xb)
        ends_s = segment_ends(el, s)
        sec = el.sections[s]

        q = fixed_end_forces(subload, sec, ends_s, xa, xb, el.rollangle)
        q̃ = condense_fef(q, sec, Ls, ends_s)

        cache.q_local[el.index][s] .+= q̃

        Qg = _rotate12(Λ', q̃)
        slots = segment_slots(el, s)
        @inbounds for i in 1:12
            cache.Pf[g[slots[i]]] += Qg[i]
        end
    end
end

# restrict a whole-member load to segment [ta, tb], re-expressed in
# segment-local fractions; nothing if the segment carries no load
function _restrict_load(load::DistributedLoad{T}, ta::Real, tb::Real, el) where {T}
    lo = max(ta, first(load.t))
    hi = min(tb, last(load.t))
    hi - lo <= eps(T) && return nothing

    # breakpoints falling inside the overlap, plus the overlap bounds
    inner = [t for t in load.t if lo < t < hi]
    tsub = vcat(lo, inner, hi)
    wsub = [_intensity_at(load, t) for t in tsub]
    tloc = (tsub .- ta) ./ (tb - ta)
    return DistributedLoad(el, tloc, wsub, load.direction;
        coords=load.coords, id=load.id, case=load.case)
end

function _restrict_load(load::PointLoad{T}, ta::Real, tb::Real, el) where {T}
    # boundary positions attach to the segment they open (last segment closes)
    inseg = (ta <= load.position < tb) || (tb == 1 && load.position == 1)
    inseg || return nothing
    τ = clamp((load.position - ta) / (tb - ta), eps(T), 1 - eps(T))
    return PointLoad(el, τ, load.value; coords=load.coords, id=load.id, case=load.case)
end

function _restrict_load(load::PointMoment{T}, ta::Real, tb::Real, el) where {T}
    inseg = (ta <= load.position < tb) || (tb == 1 && load.position == 1)
    inseg || return nothing
    τ = clamp((load.position - ta) / (tb - ta), eps(T), 1 - eps(T))
    return PointMoment(el, τ, load.value; coords=load.coords, id=load.id, case=load.case)
end

function _restrict_load(load::SelfWeight{T}, ta::Real, tb::Real, el) where {T}
    return SelfWeight(el; g=load.g, factor=load.factor, id=load.id, case=load.case)
end

# piecewise-linear intensity of a DistributedLoad at whole-member fraction t
function _intensity_at(load::DistributedLoad{T}, t::Real) where {T}
    t <= first(load.t) && return t == first(load.t) ? first(load.w) : zero(T)
    t >= last(load.t) && return t == last(load.t) ? last(load.w) : zero(T)
    s = searchsortedlast(load.t, t)
    s == length(load.t) && return last(load.w)
    frac = (t - load.t[s]) / (load.t[s+1] - load.t[s])
    return load.w[s] + (load.w[s+1] - load.w[s]) * frac
end

"""
    _rotate12(A, v) -> SVector{12}

Blockwise rotation of a 12-vector by a 3×3 matrix (four 3-blocks: start
force, start moment, end force, end moment). Written as one immutable
`SVector` expression (no scratch buffer) so it is differentiable on the pure
path; the per-entry arithmetic order is the explicit `a·x + b·y + c·z`, the
same as the in-place version it replaced, so stateful results are unchanged.
"""
function _rotate12(A::SMatrix{3,3}, v::AbstractVector)
    r(b, row) = A[row, 1] * v[3b+1] + A[row, 2] * v[3b+2] + A[row, 3] * v[3b+3]
    return SVector(
        r(0, 1), r(0, 2), r(0, 3), r(1, 1), r(1, 2), r(1, 3),
        r(2, 1), r(2, 2), r(2, 3), r(3, 1), r(3, 2), r(3, 3),
    )
end
