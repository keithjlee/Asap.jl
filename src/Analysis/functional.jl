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
    EA::Any                    # optional plain per-element axial-rigidity vector:
    # when provided, batched paths index THIS instead of mapping EA(section)
    # over the section structs — per-element struct-getfield pullbacks are
    # the single most expensive pattern under reverse AD
    ends::Any                  # optional per-element EndConditions override:
    # when provided, frame kernels take end springs from HERE instead of the
    # cache — this is what makes CONNECTION STIFFNESS a differentiable
    # design variable (semi-rigid joint optimization)
end

ModelState{T}(X::Matrix, sections::AbstractVector) where {T} =
    ModelState{T}(X, sections, nothing, nothing)
ModelState{T}(X::Matrix, sections::AbstractVector, EA) where {T} =
    ModelState{T}(X, sections, EA, nothing)

"""
    extract_state(model) -> ModelState

The model's current geometry and sections as a differentiable state.
"""
function extract_state(model::Model{T}) where {T}
    sections = Any[_state_sections(el) for el in model.elements]
    ea = [s isa AbstractVector ? zero(T) : T(EA(s)) for s in sections]   # (super-elements: unused)
    return ModelState{T}(reduce(hcat, (collect(n.position) for n in model.nodes)),
        sections, ea, nothing)
end

_state_sections(el::AbstractElement) = el.section
_state_sections(el::VariableElement) = el.sections

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
        _group_entries_all(g, state)
    end
    return reduce(vcat, parts)
end

# default: per-element kernel evaluation
function _group_entries_all(g::ElementGroup, state::ModelState)
    reduce(vcat, map(eachindex(g.elements)) do e
        _group_entries(g, e, state.sections[g.model_indices[e]],
            _position(state.X, g.i1[e]), _position(state.X, g.i2[e]),
            state.ends === nothing ? nothing : state.ends[g.model_indices[e]])
    end)
end

"""
Truss groups assemble BATCHED: one incidence matmul gives every element
vector, and the 36 stiffness entries per element are broadcasted array
expressions — no per-element closures, so reverse-mode AD pulls back
through a handful of dense array operations instead of thousands of
scalar-graph nodes. (Entry r = (b−1)·6 + a of column e is
±(EA/L)·n_{ia}·n_{ib} — exactly `vec` of the 6×6 kernel, in scatter order.)
"""
function _group_entries_all(g::ElementGroup{<:TrussElement}, state::ModelState)
    V = state.X * g.Cinc'                             # 3 × nel element vectors
    L = sqrt.(sum(abs2, V; dims=1))                   # 1 × nel
    N = V ./ L                                        # unit vectors
    ea = state.EA === nothing ? map(i -> EA(state.sections[i]), g.model_indices) :
         state.EA[g.model_indices]
    c = reshape(ea, 1, :) ./ L                        # 1 × nel  (EA/L)

    comps = (N[1:1, :], N[2:2, :], N[3:3, :])         # three row slices, once
    rows = map(1:36) do r
        a = (r - 1) % 6 + 1
        b = (r - 1) ÷ 6 + 1
        σ = (a <= 3) == (b <= 3) ? 1.0 : -1.0
        σ .* c .* comps[(a-1)%3+1] .* comps[(b-1)%3+1]
    end
    return vec(reduce(vcat, rows))
end

# Per-element active stiffness entries (scatter order), consuming ONLY plain
# data — no mutable struct reads inside the differentiated closure. Full
# frames skip slicing entirely, avoiding the expensive generic getindex
# pullback on SMatrix under Zygote.
_group_entries(::ElementGroup{<:TrussElement}, e::Int, section, x1, x2, ends) =
    vec(truss_stiffness(section, x1, x2))

function _group_entries(g::ElementGroup{<:FrameElement}, e::Int, section, x1, x2, ends)
    ec = ends === nothing ? g.endss[e]::EndConditions : ends
    Ke = frame_stiffness(section, ec, x1, x2, g.rollangles[e])
    slots = g.slots[e]
    return Base.length(slots) == 12 ? vec(Ke) : vec(Ke[slots, slots])
end

function _group_entries(g::ElementGroup{<:VariableElement}, e::Int, sections, x1, x2, ends)
    Ke = stiffness(g.elements[e], sections, x1, x2)   # super-elements keep the struct path
    return vec(Ke[g.slots[e], g.slots[e]])
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
solve_free(K::SparseMatrixCSC, F::AbstractMatrix) = _factorize(K) \ F   # multi-RHS

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
