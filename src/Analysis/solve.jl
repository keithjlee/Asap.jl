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
    offset = 6 * length(model.nodes)          # internal blocks follow nodal slots
    for (i, el) in enumerate(model.elements)
        el.index = i
        validate_rigidities(el)
        if el isa VariableElement
            el.internal_offset = offset
            offset += n_internal_dofs(el)
        end
    end
    model.cache = build_cache(model)
    return model
end

"""
    validate_rigidities(el)

Catch zero-rigidity/connection mismatches at `process!` time, before they
become an opaque singular factorization — or worse, a model that solves but
silently carries no bending in a member you meant to be a beam.

The check is exact: a rigidity may be zero ONLY if the matching end
connection is fully released at both ends (spring stiffness 0), so the
element never promises stiffness it cannot deliver:

- `EA == 0` is always an error (a member with no axial rigidity does nothing).
- `EIx == 0` requires both `kz` springs = 0 (local x–y bending released).
- `EIy == 0` requires both `ky` springs = 0 (local x–z bending released).
- `GJ == 0` requires both `kt` springs = 0 (torsion released).

The legitimate axial-only patterns pass untouched: `Section(mat, A)` on a
`TrussElement`, or on a fully released (`:freefree`) frame brace. A
`VariableElement` must have all rigidities positive in every segment —
its interior joints are rigid, so a zero-rigidity segment always creates
an unstiffened internal DOF.
"""
function validate_rigidities(el::FrameElement)
    s = el.section
    EA(s) > 0 || throw(ArgumentError(
        "element $(el.index) (:$(el.id)) has EA = 0 — a member with no axial rigidity carries nothing"))
    e1, e2 = el.ends.e1, el.ends.e2
    _plane_ok(EIx(s), e1.kz, e2.kz) || throw(ArgumentError(
        "element $(el.index) (:$(el.id)) has EIx = 0 but a moment connection in the local x–y plane " *
        "(kz ≠ 0 at an end) — use a TrussElement, release the ends, or give the section Ix"))
    _plane_ok(EIy(s), e1.ky, e2.ky) || throw(ArgumentError(
        "element $(el.index) (:$(el.id)) has EIy = 0 but a moment connection in the local x–z plane " *
        "(ky ≠ 0 at an end) — use a TrussElement, release the ends, or give the section Iy"))
    _plane_ok(GJ(s), e1.kt, e2.kt) || throw(ArgumentError(
        "element $(el.index) (:$(el.id)) has GJ = 0 but a torsional connection (kt ≠ 0 at an end) — " *
        "release torsion or give the section J"))
    return nothing
end

function validate_rigidities(el::TrussElement)
    EA(el.section) > 0 || throw(ArgumentError(
        "element $(el.index) (:$(el.id)) has EA = 0 — a member with no axial rigidity carries nothing"))
    return nothing
end

function validate_rigidities(el::VariableElement)
    for (i, s) in enumerate(el.sections)
        (EA(s) > 0 && EIx(s) > 0 && EIy(s) > 0 && GJ(s) > 0) || throw(ArgumentError(
            "VariableElement $(el.index) (:$(el.id)) segment $i has a zero rigidity — " *
            "interior joints are rigid, so every segment needs EA, EIx, EIy, GJ > 0"))
    end
    return nothing
end

validate_rigidities(::AbstractElement) = nothing

# a zero rigidity is acceptable only when the matching plane is fully released
_plane_ok(rigidity, k1, k2) = rigidity > 0 || (iszero(k1) && iszero(k2))

"""
    solve!(model; reprocess = false, solver = nothing) -> model

Run a linear static analysis: assemble the stiffness matrix and load
vectors in place, factorize (Cholesky, with LDLᵀ fallback for indefinite
spring cases), solve for the free displacements, and post-process element
end forces and support reactions into a fresh [`LinearResults`](@ref)
stored at `model.results`.

Repeated solves reuse the frozen sparsity pattern, all buffers, AND the
factorization's symbolic analysis (numeric-only refactorization). Pass
`reprocess = true` after topology changes.

`solver` selects the linear-solver backend. The default (`nothing`) is the
built-in CHOLMOD path and needs no extra packages. With LinearSolve.jl
loaded, any of its algorithms works — `solve!(model; solver = KrylovJL_CG())`
— and the choice is remembered on the model's cache for subsequent solves.
"""
function solve!(model::Model{T}; reprocess::Bool=false, solver=nothing) where {T}
    (model.cache === nothing || reprocess) && process!(model)
    cache = model.cache::AnalysisCache{T}

    assemble_K!(cache, model)
    assemble_loads!(cache, model)

    part = cache.partition
    F = cache.P[part.free] .- cache.Pf[part.free]

    fc = _ensure_factorization!(cache, solver)
    uf = _backsolve(fc.solver, fc, F)

    u = zeros(T, part.n_global)
    u[part.free] = uf

    model.results = _post_process(cache, model, u, dot(uf, F))
    return model
end

"""
    FactorizationCache

Holder for the factorization state behind the solver seam: which `solver`
backend produced it (`nothing` = the built-in CHOLMOD path), the backend
factorization object `F`, and backend-private `meta` (the built-in path
stores whether it fell back to LDLᵀ). Opaque solver state — declared
tangent-free for AD engines, so ANY backend's internals stay invisible to
differentiation.
"""
mutable struct FactorizationCache
    solver::Any
    F::Any
    meta::Any
end

#=
The solver seam. Three functions, dispatched on the solver object:

    _factorize(solver, K)      -> FactorizationCache   (symbolic + numeric)
    _refactorize!(solver, fc, K) -> fc                  (numeric only — the
                                    sparsity pattern is frozen, so repeated
                                    solves reuse the symbolic analysis)
    _backsolve(solver, fc, b)  -> x                     (vector or multi-RHS)

`solver === nothing` is the zero-dependency default (CHOLMOD Cholesky with
LDLᵀ fallback). Extensions add methods for their own solver types — e.g.
`AsapLinearSolveExt` accepts any `LinearSolve.SciMLLinearSolveAlgorithm`:

    using LinearSolve
    solve!(model; solver = KrylovJL_CG())

The chosen solver is remembered on the cache: subsequent `solve!` calls
reuse it (and its symbolic analysis) without re-passing the keyword.
=#
function _factorize(::Nothing, K::SparseMatrixCSC)
    try
        return FactorizationCache(nothing, cholesky(Symmetric(K)), false)
    catch err
        err isa LinearAlgebra.PosDefException || rethrow()
        return FactorizationCache(nothing, ldlt(Symmetric(K)), true)
    end
end

function _refactorize!(::Nothing, fc::FactorizationCache, K::SparseMatrixCSC)
    try
        if fc.meta === true
            ldlt!(fc.F, Symmetric(K))
        else
            cholesky!(fc.F, Symmetric(K))
        end
    catch err
        err isa LinearAlgebra.PosDefException || rethrow()
        newfc = _factorize(nothing, K)
        fc.F = newfc.F
        fc.meta = newfc.meta
    end
    return fc
end

_backsolve(::Nothing, fc::FactorizationCache, b::AbstractVecOrMat) = fc.F \ b

# convenience: the wrapper solves like a factorization (`fc \ b`)
Base.:\(fc::FactorizationCache, b::AbstractVecOrMat) = _backsolve(fc.solver, fc, b)

_factorize(solver, K::SparseMatrixCSC) =
    error("no solver backend for $(typeof(solver)) — load LinearSolve.jl " *
          "(or define the _factorize/_refactorize!/_backsolve seam methods for it)")

"""
    _ensure_factorization!(cache, solver) -> FactorizationCache

Reuse the cache's factorization when possible: same (or unspecified)
solver → numeric-only refactorization on the frozen pattern; first solve or
a different solver → fresh symbolic + numeric factorization, remembered for
subsequent solves.
"""
function _ensure_factorization!(cache::AnalysisCache, solver)
    fc = cache.factorization
    if fc isa FactorizationCache && (solver === nothing || solver === fc.solver)
        return _refactorize!(fc.solver, fc, cache.K)
    end
    fc = _factorize(solver, cache.K)
    cache.factorization = fc
    return fc
end

# raw single-shot factorization (pure path / rrule use — no cache, default
# backend): Cholesky when SPD, LDLᵀ otherwise
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
        forces[el.index] = _recover_forces!(reactions, isfixed, cache, el, u)
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

# primitive elements: one 12-vector of local end forces
function _recover_forces!(reactions, isfixed, cache::AnalysisCache{T},
    el::Union{FrameElement,TrussElement}, u::Vector{T}) where {T}
    x1 = el.nodeStart.position
    x2 = el.nodeEnd.position
    Ke = stiffness(el, x1, x2)
    g = element_global_dofs(el)
    u_el = SVector{12,T}(ntuple(i -> u[g[i]], 12))

    Λ = local_frame(x1, x2, el isa FrameElement ? el.rollangle : zero(T))
    Qg = _rotate12(Λ', SVector{12,T}(cache.q_local[el.index][1]))

    f_global = Ke * u_el + Qg
    @inbounds for i in 1:12
        isfixed[g[i]] && (reactions[g[i]] += f_global[i])
    end
    return Vector{T}(_rotate12(Λ, f_global))
end

# VariableElement: per-segment recovery (each segment behaves like a
# primitive element over its slice of the element's DOF layout); returns
# the concatenated per-segment local end-force vectors (12 per segment)
function _recover_forces!(reactions, isfixed, cache::AnalysisCache{T},
    el::VariableElement, u::Vector{T}) where {T}
    g = element_global_dofs(el)
    Λ = local_frame(el)
    out = zeros(T, 12 * n_segments(el))

    for s in 1:n_segments(el)
        xa, xb = segment_positions(el, s)
        Ks = frame_stiffness(el.sections[s], segment_ends(el, s), xa, xb, el.rollangle)
        slots = segment_slots(el, s)
        u_s = SVector{12,T}(ntuple(i -> u[g[slots[i]]], 12))
        Qg = _rotate12(Λ', SVector{12,T}(cache.q_local[el.index][s]))

        f_global = Ks * u_s + Qg
        @inbounds for i in 1:12
            gi = g[slots[i]]
            isfixed[gi] && (reactions[gi] += f_global[i])
        end
        out[12*(s-1).+(1:12)] = Vector{T}(_rotate12(Λ, f_global))
    end
    return out
end
