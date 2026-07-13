#=
The analysis pipeline: process! (once per topology) → solve! (per state).
=#

"""
    process!(model) -> model

Build the model's analysis structure: assign node/element indices, classify
DOFs (free / fixed / inactive), group elements by type, and freeze the
free×free sparsity pattern with its scatter maps (see [`AnalysisCache`](@ref)).

Call once per *topology*. Geometry, section, spring-value, and load changes
do NOT require re-processing — they are picked up by the next `solve!`.
Changes that alter connectivity or DOF activity (adding/removing
nodes/elements/springs, changing end conditions between released and
engaged) do.
"""
function process!(model::Model{T}) where {T}
    for (i, n) in enumerate(model.nodes)
        n.index = i
    end
    for (i, el) in enumerate(model.elements)
        el.index = i
    end
    model.cache = build_cache(model)
    return model
end

"""
    solve!(model; reprocess = false) -> model

Run a linear static analysis: assemble the stiffness matrix and load
vectors in place, factorize (Cholesky, with LDLᵀ fallback for indefinite
spring cases), solve for the free displacements, and post-process element
end forces and support reactions into a fresh [`LinearResults`](@ref)
stored at `model.results`.

Repeated solves reuse the frozen sparsity pattern and all buffers. Pass
`reprocess = true` after topology changes.
"""
function solve!(model::Model{T}; reprocess::Bool=false) where {T}
    (model.cache === nothing || reprocess) && process!(model)
    cache = model.cache::AnalysisCache{T}

    assemble_K!(cache, model)
    assemble_loads!(cache, model)

    part = cache.partition
    F = cache.P[part.free] .- cache.Pf[part.free]

    cache.factorization = _factorize(cache.K)
    uf = cache.factorization \ F

    u = zeros(T, part.n_global)
    u[part.free] = uf

    model.results = _post_process(cache, model, u, dot(uf, F))
    return model
end

# factorization seam: Cholesky when SPD, LDLᵀ otherwise (kept behind one
# function so a LinearSolve.jl extension can swap strategies later)
function _factorize(K::SparseMatrixCSC)
    try
        return cholesky(Symmetric(K))
    catch err
        err isa LinearAlgebra.PosDefException || rethrow()
        return ldlt(Symmetric(K))
    end
end

"""
Recover element end forces and support reactions element-wise — no
full-space stiffness matrix is ever formed. For each element:
`f_global = Kₑ·u_el + Q_global`; the local end-force vector is its
blockwise rotation into element coordinates, and its fixed-DOF entries
accumulate into the reactions. Spring reactions (−k·u) accumulate at
supported spring DOFs.
"""
function _post_process(cache::AnalysisCache{T}, model::Model{T}, u::Vector{T},
    compliance::T) where {T}

    reactions = zeros(T, cache.partition.n_global)
    isfixed = falses(cache.partition.n_global)
    isfixed[cache.partition.fixed] .= true

    forces = Vector{Vector{T}}(undef, length(model.elements))
    for el in model.elements
        x1 = el.nodeStart.position
        x2 = el.nodeEnd.position
        Ke = stiffness(el, x1, x2)
        g = element_global_dofs(el)
        u_el = SVector{12,T}(u[g[1]], u[g[2]], u[g[3]], u[g[4]], u[g[5]], u[g[6]],
            u[g[7]], u[g[8]], u[g[9]], u[g[10]], u[g[11]], u[g[12]])

        q̃ = cache.q_local[el.index]
        Λ = local_frame(x1, x2, el isa FrameElement ? el.Ψ : zero(T))
        Qg = _rotate12(Λ', SVector{12,T}(q̃))

        f_global = Ke * u_el + Qg
        forces[el.index] = Vector{T}(_rotate12(Λ, f_global))

        @inbounds for i in 1:12
            isfixed[g[i]] && (reactions[g[i]] += f_global[i])
        end
    end

    for sp in model.springs
        s = node_dof_start(sp.node)
        @inbounds for i in 1:6
            if sp.stiffness[i] > 0 && isfixed[s+i]
                reactions[s+i] -= sp.stiffness[i] * u[s+i]
            end
        end
    end

    return LinearResults{T}(u, reactions, forces, compliance)
end
