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
    @inbounds for (pos, k) in cache.spring_nz
        nz[pos] += k
    end
    return cache.K
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
function assemble_loads!(cache::AnalysisCache{T}, model::Model{T}) where {T}
    fill!(cache.P, zero(T))
    fill!(cache.Pf, zero(T))
    for q in cache.q_local
        fill!(q, zero(T))
    end

    inactive = Set(cache.partition.inactive)

    for load in model.loads
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
    el isa FrameElement || error("element loads require a FrameElement " *
                                 "(a $(nameof(typeof(el))) cannot carry transverse load) — " *
                                 "apply NodeForces or subdivide instead")
    x1 = el.nodeStart.position
    x2 = el.nodeEnd.position
    L = element_length(x1, x2)

    q = fixed_end_forces(load, el.section, el.ends, x1, x2, el.Ψ)
    q̃ = condense_fef(q, el.section, L, el.ends)

    cache.q_local[el.index] .+= q̃

    Λ = local_frame(x1, x2, el.Ψ)
    Qg = _rotate12(Λ', q̃)                     # local → global, blockwise
    g = element_global_dofs(el)
    @inbounds for i in 1:12
        cache.Pf[g[i]] += Qg[i]
    end
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
