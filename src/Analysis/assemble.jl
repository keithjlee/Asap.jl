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
"""
function assemble_loads!(cache::AnalysisCache{T}, model::Model{T};
    case::Union{Nothing,Symbol}=nothing) where {T}
    fill!(cache.P, zero(T))
    fill!(cache.Pf, zero(T))
    for qs in cache.q_local, q in qs
        fill!(q, zero(T))
    end

    inactive = Set(cache.partition.inactive)

    for load in model.loads
        case === nothing || load.case === case || continue
        _apply_load!(cache, load, inactive)
    end
    return cache.P, cache.Pf
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
    L = element_length(x1, x2)

    q = fixed_end_forces(load, el.section, el.ends, x1, x2, el.rollangle)
    q̃ = condense_fef(q, el.section, L, el.ends)

    cache.q_local[el.index][1] .+= q̃

    Λ = local_frame(x1, x2, el.rollangle)
    Qg = _rotate12(Λ', q̃)                     # local → global, blockwise
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

# blockwise rotation of a 12-vector by a 3×3 matrix (four 3-blocks)
function _rotate12(A::SMatrix{3,3}, v::AbstractVector)
    T = promote_type(eltype(A), eltype(v))
    out = _mutable12vec(T)
    @inbounds for b in 0:3
        i = 3b
        out[i+1] = A[1, 1] * v[i+1] + A[1, 2] * v[i+2] + A[1, 3] * v[i+3]
        out[i+2] = A[2, 1] * v[i+1] + A[2, 2] * v[i+2] + A[2, 3] * v[i+3]
        out[i+3] = A[3, 1] * v[i+1] + A[3, 2] * v[i+2] + A[3, 3] * v[i+3]
    end
    return SVector{12,T}(out)
end
