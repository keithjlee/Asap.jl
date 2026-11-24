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

Loads: nodal loads `P` are constant data read from the cache. Element-load
fixed-end forces are DESIGN-DEPENDENT: the pure path re-evaluates them at the
state's positions, sections and end conditions ([`assemble_Pf`](@ref)), so
gradients include the load path (a member's line load scales with its
length, self-weight with ρA, condensed end moments with the end springs and
EI). Pass `loads = :frozen` to [`solve`](@ref)/[`compliance`](@ref) for the
cached, design-constant fixed-end forces instead.
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
    #pure-path spring diagonals read the model's CURRENT springs (same
    #freshness contract as the in-place assembler)
    model.cache === nothing || _refresh_spring_nzvec!(model.cache, model)
    sections = Any[_state_sections(el) for el in model.elements]
    ea = [s isa AbstractSection ? T(EA(s)) : zero(T) for s in sections]  # (super-elements, shells: unused)
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

# generic N-node default (surface elements and other non-line types):
# per-element kernel evaluation at the state's positions of ALL its nodes.
# Element types specialize `_group_entries_nnode` to consume plain data only
# (the default reads the element struct — fine for prototypes, not for
# reverse-mode AD at scale).
function _group_entries_all(g::ElementGroup, state::ModelState)
    return _group_entries_nnode_all(g, state, Val(size(g.inodes, 1)))
end

function _group_entries_nnode_all(g::ElementGroup, state::ModelState, ::Val{N}) where {N}
    reduce(vcat, map(eachindex(g.elements)) do e
        X = ntuple(k -> _position(state.X, g.inodes[k, e]), Val(N))
        _group_entries_nnode(g, e, state.sections[g.model_indices[e]], X)
    end)
end

function _group_entries_nnode(g::ElementGroup, e::Int, section, X::NTuple)
    Ke = stiffness(g.elements[e], X...)
    slots = g.slots[e]
    return Base.length(slots) == size(Ke, 1) ? vec(Ke) : vec(Ke[slots, slots])
end

# shells: the kernel takes only (section, positions) — no element-struct reads
# inside the differentiated closure; every slot is active (no slicing)
_group_entries_nnode(::ElementGroup{<:ShellElement}, e::Int, section, X::NTuple) =
    vec(shell_stiffness(section, X...))

# line elements: per-element kernel evaluation from the two endpoint columns
function _group_entries_all(g::ElementGroup{<:AbstractLineElement}, state::ModelState)
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
    solve_free(solver, K, F) -> u_free

Solve the reduced linear system. The AD extension supplies the adjoint
(one extra back-substitution on the cached factorization — the classical
adjoint method).

The 3-arg form goes through the solver seam (`_factorize(solver, K)`):
`solver` may be `nothing` (built-in CHOLMOD), any LinearSolve algorithm
(with `AsapLinearSolveExt`), or a [`CachedSolver`](@ref) wrapper that
reuses one factorization across value/adjoint/tangent passes at the same
design. The 2-arg form is the default-solver fast path — kept as its own
method (not a forwarder) because Enzyme's imported rules match this exact
signature.
"""
solve_free(K::SparseMatrixCSC, F::AbstractVector) = _factorize(K) \ F
solve_free(K::SparseMatrixCSC, F::AbstractMatrix) = _factorize(K) \ F   # multi-RHS
solve_free(solver, K::SparseMatrixCSC, F::AbstractVecOrMat) = _factorize(solver, K) \ F

"""
    assemble_Pf(cache, state) -> Vector

The full-space fixed-end force vector `Pf` evaluated at the design `state`:
every element load applied by the last [`assemble_loads!`](@ref) (the
plain-data [`FEFTable`](@ref) in `cache.fef`) is lowered through
[`element_fef_global`](@ref) at the state's node positions, the element's
state section, and its end conditions (`state.ends` when the state carries
them, else the model's), then scattered to global DOFs.

Pure and differentiable: one map over the loaded elements of immutable
12-vector kernels and one sparse matvec. With no element loads it returns
the cached `Pf` untouched, so nodal-load-only models pay nothing. Loads on
`VariableElement`s enter as the constant `cache.fef.static` (their internal
geometry is not part of the state).

At the model's own state (`extract_state(model)`) this reproduces the
stateful `cache.Pf` to rounding.
"""
function assemble_Pf(cache::AnalysisCache, state::ModelState)
    t = cache.fef
    isempty(t.loads) && isempty(t.sloads) && return cache.Pf
    isempty(t.sloads) && return t.static + _frame_fef(t, state)
    isempty(t.loads) && return t.static + _shell_fef(t, state)
    return t.static + _frame_fef(t, state) + _shell_fef(t, state)
end

function _frame_fef(t::FEFTable, state::ModelState)
    Q = map(eachindex(t.loads)) do j
        e = t.elem[j]
        ends = state.ends === nothing ? t.ends[j] : state.ends[e]
        element_fef_global(t.loads[j], state.sections[e], ends,
            _position(state.X, t.i1[j]), _position(state.X, t.i2[j]), t.rollangle[j])
    end
    return t.scatter * reduce(vcat, Q)
end

# shell surface loads at the state: traction (local pressure follows the
# element orientation; self-weight follows ρt of the state section) →
# condensed consistent load → global 24-vector, scattered
function _shell_fef(t::FEFTable, state::ModelState)
    Q = map(eachindex(t.sloads)) do j
        sec = state.sections[t.selem[j]]
        X = (_position(state.X, t.snodes[1, j]), _position(state.X, t.snodes[2, j]),
             _position(state.X, t.snodes[3, j]), _position(state.X, t.snodes[4, j]))
        _shell_fef_global(t.sloads[j], sec, X)
    end
    return t.sscatter * reduce(vcat, Q)
end

function _shell_fef_global(load, sec, X)
    R = shell_frame(X...)[1]
    q = _shell_traction_state(load, sec, R)
    return _shell_fef_plain(membrane_rigidity(sec), q, X...)
end

# global fixed-end 24-vector of a uniform global traction q: plain-data AD target
function _shell_fef_plain(A, q, X1, X2, X3, X4)
    R, xy, _ = shell_frame(X1, X2, X3, X4)
    f24, fb = shell_consistent_load(xy, R * q)
    K16 = _itw_K16(xy, A)
    b = SVector(13, 14, 15, 16)
    corr = K16[_SONE12, b] * (K16[b, b] \ fb)
    Δ = vcat(corr, zero(eltype(corr)))[_IDX_MEMBRANE]
    return _rotate24(R', -(f24 - Δ))
end

"""
    fef_local(cache, state) -> Matrix   (12 × n_loads)

The LOCAL condensed fixed-end force vectors `q̃` ([`condensed_fef`](@ref)) of
every state-tracked element load at the design `state`, one column per load,
in `cache.fef` order (`cache.fef.elem[j]` is load `j`'s element). These are
the `q̃` terms of the local force recovery `f = k·u_local + q̃`, e.g. slot 7
is the load's contribution to the end-node axial action. Pure and
differentiable.
"""
function fef_local(cache::AnalysisCache, state::ModelState)
    t = cache.fef
    T = eltype(state.X)
    isempty(t.loads) && return zeros(T, 12, 0)
    Q = map(eachindex(t.loads)) do j
        e = t.elem[j]
        ends = state.ends === nothing ? t.ends[j] : state.ends[e]
        condensed_fef(t.loads[j], state.sections[e], ends,
            _position(state.X, t.i1[j]), _position(state.X, t.i2[j]), t.rollangle[j])
    end
    return reshape(reduce(vcat, Q), 12, :)
end

# does the pure path re-evaluate fixed-end forces? (validates `loads`)
function _design_loads(cache::AnalysisCache, loads::Symbol)
    loads === :design || loads === :frozen ||
        throw(ArgumentError("loads must be :design or :frozen, got :$loads"))
    return loads === :design && !(isempty(cache.fef.loads) && isempty(cache.fef.sloads))
end

# the reduced solve, behind a function barrier: each caller branch passes a
# concretely-typed F, so the solve_free call site stays type-stable on the
# constant-load path (Enzyme's imported rules and forward mode need that;
# the design-load path maps over abstract load objects and may infer loosely)
function _pure_solve(K, F, solver)
    # the default keeps the 2-arg call — Enzyme's imported rules match it
    return solver === nothing ? solve_free(K, F) : solve_free(solver, K, F)
end

# reduced solve through the diaphragm map: u_free = C·(CᵀKC)⁻¹·CᵀF — K̂'s values
# are a sparse linear map of K's (differentiable); C is frozen at the last
# solve!/process! (offsets constant under AD)
function _constrained_solve(cm::ConstraintMap, K, F, solver)
    Khat = _reduce_K(cm, K)
    return cm.C * _pure_solve(Khat, cm.C' * F, solver)
end

"""
    solve(model, state; solver = nothing, loads = :design) -> Vector

Pure linear solve on a processed model: assemble from the differentiable
`state`, solve, and return the FULL-space displacement vector (inactive and
fixed slots zero). Requires `process!` (or a prior `solve!`) to have built
the cache; nodal loads are taken from the cache as constants, so call
`assemble_loads!` first if loads changed.

`loads = :design` (default) re-evaluates element-load fixed-end forces at
the state ([`assemble_Pf`](@ref)): a member's line load follows its current
length, self-weight its ρA, condensed end moments its end springs and EI.
`loads = :frozen` uses the cached, design-constant `Pf` instead (the
behavior before v1.2: loads fixed at the model's stored geometry).

`solver` selects the linear-solver backend for the reduced solve (see
[`solve_free`](@ref)); `nothing` is the built-in CHOLMOD path.

This is the AD entry point: gradients of any scalar of the returned vector
with respect to node positions and section properties flow through Zygote
with no further ceremony.
"""
function solve(model::Model{T}, state::ModelState; solver = nothing, loads::Symbol = :design) where {T}
    cache = model.cache::AnalysisCache{T}
    K = assemble_K(cache, state)
    free = cache.partition.free
    if cache.constraint !== nothing
        F = _design_loads(cache, loads) ? (cache.P - assemble_Pf(cache, state))[free] : (cache.P - cache.Pf)[free]
        return cache.free_embed * _constrained_solve(cache.constraint, K, F, solver)
    end
    if _design_loads(cache, loads)
        uf = _pure_solve(K, (cache.P - assemble_Pf(cache, state))[free], solver)
    else
        uf = _pure_solve(K, (cache.P - cache.Pf)[free], solver)
    end
    return cache.free_embed * uf
end

"""
    compliance(model, state; solver = nothing, loads = :design) -> T

External work `Fᵀu_free` of the pure solve — the canonical smooth stiffness
objective, differentiable end-to-end. `loads` as in [`solve`](@ref): with
design-dependent loads `F` itself varies with the design, and the gradient
includes that term.
"""
function compliance(model::Model{T}, state::ModelState; solver = nothing, loads::Symbol = :design) where {T}
    cache = model.cache::AnalysisCache{T}
    K = assemble_K(cache, state)
    free = cache.partition.free
    if cache.constraint !== nothing
        F = _design_loads(cache, loads) ? (cache.P - assemble_Pf(cache, state))[free] : (cache.P - cache.Pf)[free]
        return sum(F .* _constrained_solve(cache.constraint, K, F, solver))
    end
    if _design_loads(cache, loads)
        F = (cache.P - assemble_Pf(cache, state))[free]
        # sum∘broadcast, not LinearAlgebra.dot: dot lowers to BLAS, which
        # Enzyme's forward mode cannot handle under runtime activity
        return sum(F .* _pure_solve(K, F, solver))
    else
        F0 = (cache.P - cache.Pf)[free]
        return sum(F0 .* _pure_solve(K, F0, solver))
    end
end
