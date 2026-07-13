"""
    ModelState{T}

The differentiable state of a model: everything automatic differentiation
is allowed to vary, extracted as plain data. Element kernels are pure
functions, so evaluating the pipeline on a `ModelState` involves **no
mutation anywhere** — Zygote (and other AD engines) differentiate it
directly, with only two custom rules (sparse construction and the linear
solve) supplied by the ChainRules extension.

# Fields
- `X::Matrix{T}`: node positions as a `3 × n_nodes` matrix (column `i` is
  node `i`) — plain arrays differentiate cleanly through every AD engine;
  static vectors remain internal to the kernels
- `sections::Vector{AbstractSection{T}}`: per-ELEMENT sections, indexed
  like `model.elements` (elements sharing a section object simply repeat it)

Construct the reference state with [`extract_state`](@ref), then produce
perturbed states in your objective (e.g. rebuild positions from a design
vector, or sections from area variables) — see the gradient tests for
patterns.

Loads are treated as constant data in the pure path for now (`P` and `Pf`
are read from the cache); differentiating through load values/FEFs is a
Phase 5a extension if needed.
"""
struct ModelState{T}
    X::Matrix{T}
    sections::AbstractVector   # per-element AbstractSections; loosely typed so
    # AD-constructed section vectors (e.g. from map over a design vector)
    # pass through without conversion (conversion mutates — Zygote-hostile)
end

"""
    extract_state(model) -> ModelState

The model's current geometry and sections as a differentiable state.
"""
extract_state(model::Model{T}) where {T} =
    ModelState{T}(reduce(hcat, (collect(n.position) for n in model.nodes)),
        AbstractSection{T}[el.section for el in model.elements])

# column i of the position matrix as a static vector (pure expression —
# Zygote differentiates the scalar getindex calls cleanly)
@inline _position(X::AbstractMatrix, i::Int) = SVector(X[1, i], X[2, i], X[3, i])

"""
    stiffness_entries(cache, state) -> Vector

All element stiffness entries in scatter order (group → element →
column-major active-slot pairs): the differentiable vector `V` with
`nzval(K) = cache.scatter · V + cache.spring_nzvec`. Pure — safe under
Zygote.
"""
function stiffness_entries(cache::AnalysisCache, state::ModelState)
    parts = map(cache.groups) do g
        reduce(vcat, map(eachindex(g.elements)) do e
            el = g.elements[e]
            Ke = stiffness(el, state.sections[g.model_indices[e]],
                _position(state.X, el.nodeStart.index), _position(state.X, el.nodeEnd.index))
            vec(Ke[g.slots[e], g.slots[e]])
        end)
    end
    return reduce(vcat, parts)
end

"""
    assemble_K(cache, state) -> SparseMatrixCSC

Pure (non-mutating) assembly of the free×free stiffness matrix from a
differentiable [`ModelState`](@ref). Identical numerics to
[`assemble_K!`](@ref) — the two paths share every element kernel and the
frozen sparsity pattern — verified by parity tests.
"""
function assemble_K(cache::AnalysisCache, state::ModelState)
    V = stiffness_entries(cache, state)
    nz = cache.scatter * V + cache.spring_nzvec
    return sparse_from_pattern(cache.K, nz)
end

"""
    sparse_from_pattern(pattern, nzval) -> SparseMatrixCSC

A sparse matrix with `pattern`'s frozen structure and the given values.
The AD extension supplies its reverse rule (cotangents project onto the
pattern's stored entries).
"""
sparse_from_pattern(pattern::SparseMatrixCSC, nzval::AbstractVector) =
    SparseMatrixCSC(size(pattern, 1), size(pattern, 2),
        copy(pattern.colptr), copy(pattern.rowval), nzval)

"""
    solve_free(K, F) -> u_free

Solve the reduced linear system. The AD extension supplies the adjoint
(one extra back-substitution on the cached factorization — the classical
adjoint method).
"""
solve_free(K::SparseMatrixCSC, F::AbstractVector) = _factorize(K) \ F

"""
    solve(model, state; loads = :cached) -> Vector

Pure linear solve on a processed model: assemble from the differentiable
`state`, solve, and return the FULL-space displacement vector (inactive and
fixed slots zero). Requires `process!` (or a prior `solve!`) to have built
the cache; load vectors are taken from the cache as constants — call
`assemble_loads!` first if loads changed.

This is the AD entry point: gradients of any scalar of the returned vector
with respect to node positions and section properties flow through Zygote
with no further ceremony.
"""
function solve(model::Model{T}, state::ModelState) where {T}
    cache = model.cache::AnalysisCache{T}
    K = assemble_K(cache, state)
    F = (cache.P-cache.Pf)[cache.partition.free]
    uf = solve_free(K, F)
    return cache.free_embed * uf
end

"""
    compliance(model, state) -> T

External work `Fᵀu_free` of the pure solve — the canonical smooth stiffness
objective, differentiable end-to-end.
"""
function compliance(model::Model{T}, state::ModelState) where {T}
    cache = model.cache::AnalysisCache{T}
    K = assemble_K(cache, state)
    F = (cache.P-cache.Pf)[cache.partition.free]
    uf = solve_free(K, F)
    return dot(F, uf)
end
