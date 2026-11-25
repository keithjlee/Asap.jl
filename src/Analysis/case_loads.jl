#=
Load cases on the pure (differentiable) path.

`solve_cases!` (cases.jl) handles several cases on the stateful path. The pure
path (`solve(model, state)`) sees ONE load vector: the cache's, from the last
`assemble_loads!`. `CaseLoads` is the plain-data table that lets the pure path
solve several cases at one design with one factorization, every element load
re-evaluated at the design (design-dependent fixed-end forces), and each case
kept in its own column.

The element loads of all listed cases sit in one `FEFTable` (the same plain
data `assemble_Pf` consumes), so the per-design work is one map per load kind
over the stacked kernels. Each load belongs to exactly one case: its scatter
rows are shifted into its case's block of an `n_global·n_cases`-row CASE
SCATTER, so the fixed-end forces of every case come out of one sparse matvec
(`reshape(cscatter · q, n_global, n_cases)`), with no dense masks and no loop
over cases.
=#

"""
    CaseLoads{T}

Several load cases of a processed model as plain data, for the pure
(differentiable) path: built once by [`case_loads`](@ref), consumed by
[`solve`](@ref)`(model, state, cl)` and [`compliance`](@ref)`(model, state, cl)`
at any design state.

Column `k` of every matrix belongs to case `names[k]`; `n` is the number of
global DOFs (`6 × n_nodes`).

# Fields
- `names::Vector{Symbol}`: the case tags, in column order
- `P::Matrix{T}`: nodal loads per case (`n × n_cases`) [force, moment]
- `static::Matrix{T}`: fixed-end forces that do not follow the design (element
  loads on `VariableElement`s) per case [force, moment]
- `F0::Matrix{T}`: the frozen load vectors `P − Pf` per case, with every
  fixed-end force evaluated at the model's own geometry and sections (what
  `loads = :frozen` uses) [force, moment]
- `fef::FEFTable{T}`: the element loads of ALL listed cases (frames, quad and
  triangle shells), in case order: the plain data the design-dependent
  fixed-end forces are re-evaluated from
- `fcase`, `scase`, `tcase::Vector{Int}`: the case (column) of each frame,
  quad-shell and triangle-shell load in `fef`
- `cscatter`, `csscatter`, `ctscatter::SparseMatrixCSC{T,Int}`: CASE SCATTERS
  (`n·n_cases × 12·n_frame`, `× 24·n_quad`, `× 18·n_tri`): the ordinary
  scatters of `fef` with each load's rows moved into its case's block, so
  `reshape(cscatter * q, n, n_cases)` places every stacked fixed-end vector
  `q` in its own case's column

Built from the cache's DOF numbering: rebuild it after `process!` changes the
model's topology or loads.
"""
struct CaseLoads{T}
    names::Vector{Symbol}
    P::Matrix{T}
    static::Matrix{T}
    F0::Matrix{T}
    fef::FEFTable{T}
    fcase::Vector{Int}
    scase::Vector{Int}
    tcase::Vector{Int}
    cscatter::SparseMatrixCSC{T,Int}
    csscatter::SparseMatrixCSC{T,Int}
    ctscatter::SparseMatrixCSC{T,Int}
end

Base.length(cl::CaseLoads) = length(cl.names)

"""
    case_loads(model; cases = load_cases(model)) -> CaseLoads

The load cases `cases` of a processed model (`process!` or `solve!` first) as
a [`CaseLoads`](@ref) table for the pure path. Each case's loads are
assembled once (the same lowering as [`assemble_loads!`](@ref)`(…; case)`):
nodal loads, frozen fixed-end forces at the model's own state, and the
element loads that the pure path re-evaluates per design.

The cache is left EXACTLY as it was (its load vectors, recovery terms and
fixed-end-force table are restored), so the single-case pure path and the
stateful results are unaffected.

```julia
process!(model)
cl = case_loads(model)                        # every case, first-appearance order
cl = case_loads(model; cases = [:dead, :live])
U  = solve(model, extract_state(model), cl)   # n_dof × 2, one factorization
```
"""
function case_loads(model::Model{T}; cases::AbstractVector{Symbol} = load_cases(model)) where {T}
    cache = model.cache
    cache === nothing && throw(ArgumentError("case_loads needs a processed model: call process!(model) first"))
    cache = cache::AnalysisCache{T}
    isempty(cases) && throw(ArgumentError("case_loads needs at least one case"))
    allunique(cases) || throw(ArgumentError("case names must be distinct: $cases"))
    present = load_cases(model)
    for c in cases
        c in present || throw(ArgumentError("no load in the model carries case :$c"))
    end
    n = cache.partition.n_global
    nc = length(cases)
    P = zeros(T, n, nc)
    S = zeros(T, n, nc)
    F0 = zeros(T, n, nc)
    floads, sloads, tloads = Any[], Any[], Any[]
    fcase, scase, tcase = Int[], Int[], Int[]
    # snapshot the cache's load state; restored exactly afterwards
    P_, Pf_, fef_ = copy(cache.P), copy(cache.Pf), cache.fef
    q_ = [map(copy, qs) for qs in cache.q_local]
    try
        for (k, c) in enumerate(cases)
            assemble_loads!(cache, model; case = c)
            P[:, k] = cache.P
            S[:, k] = cache.fef.static
            F0[:, k] = cache.P .- cache.Pf
            t = cache.fef
            append!(floads, t.loads); append!(fcase, fill(k, length(t.loads)))
            append!(sloads, t.sloads); append!(scase, fill(k, length(t.sloads)))
            append!(tloads, t.tloads); append!(tcase, fill(k, length(t.tloads)))
        end
    finally
        cache.P .= P_
        cache.Pf .= Pf_
        cache.fef = fef_
        for (qs, qs_) in zip(cache.q_local, q_), (q, q0) in zip(qs, qs_)
            q .= q0
        end
    end
    fef = _build_fef_table(T, floads, n, zeros(T, n), sloads, tloads)
    return CaseLoads{T}(collect(cases), P, S, F0, fef, fcase, scase, tcase,
        _case_scatter(fef.scatter, fcase, 12, n, nc), _case_scatter(fef.sscatter, scase, 24, n, nc),
        _case_scatter(fef.tscatter, tcase, 18, n, nc))
end

# the scatter with load j's rows (columns m(j−1)+1 … mj) shifted into its case's block
function _case_scatter(scatter::SparseMatrixCSC{T,Int}, case::Vector{Int}, m::Int, n::Int, nc::Int) where {T}
    I, J, V = findnz(scatter)
    I = [i + n * (case[(j - 1) ÷ m + 1] - 1) for (i, j) in zip(I, J)]
    return sparse(I, J, V, n * nc, size(scatter, 2))
end

# the stacked (unscattered) fixed-end vectors of every load of one kind at the
# state: the maps that `_frame_fef`, `_shell_fef` and `_shell3_fef` scatter
function _frame_fef_q(t::FEFTable, state::ModelState)
    Q = map(eachindex(t.loads)) do j
        e = t.elem[j]
        ends = state.ends === nothing ? t.ends[j] : state.ends[e]
        element_fef_global(t.loads[j], state.sections[e], ends,
            _position(state.X, t.i1[j]), _position(state.X, t.i2[j]), t.rollangle[j])
    end
    return reduce(vcat, Q)
end

function _shell_fef_q(t::FEFTable, state::ModelState)
    Q = map(eachindex(t.sloads)) do j
        sec = state.sections[t.selem[j]]
        X = (_position(state.X, t.snodes[1, j]), _position(state.X, t.snodes[2, j]),
             _position(state.X, t.snodes[3, j]), _position(state.X, t.snodes[4, j]))
        _shell_fef_global(t.sloads[j], sec, X)
    end
    return reduce(vcat, Q)
end

function _shell3_fef_q(t::FEFTable, state::ModelState)
    Q = map(eachindex(t.tloads)) do j
        sec = state.sections[t.telem[j]]
        X = (_position(state.X, t.tnodes[1, j]), _position(state.X, t.tnodes[2, j]),
             _position(state.X, t.tnodes[3, j]))
        _shell_fef_global(t.tloads[j], sec, X)
    end
    return reduce(vcat, Q)
end

"""
    case_Pf(cl::CaseLoads, state) -> Matrix

The fixed-end forces of every case at the design `state` (`n × n_cases`): the
design-independent part `cl.static` plus every element load of each case
re-evaluated at the state (frame loads, quad and triangle shell surface
loads), each in its case's column. The multi-case form of
[`assemble_Pf`](@ref); pure and differentiable.
"""
function case_Pf(cl::CaseLoads, state::ModelState)
    t = cl.fef
    n, nc = size(cl.P)
    Pf = cl.static
    isempty(t.loads) || (Pf = Pf + reshape(cl.cscatter * _frame_fef_q(t, state), n, nc))
    isempty(t.sloads) || (Pf = Pf + reshape(cl.csscatter * _shell_fef_q(t, state), n, nc))
    isempty(t.tloads) || (Pf = Pf + reshape(cl.ctscatter * _shell3_fef_q(t, state), n, nc))
    return Pf
end

# per-case load vectors (full space, n × n_cases) at the state
function _case_F(cl::CaseLoads, state::ModelState, loads::Symbol)
    loads === :design || loads === :frozen ||
        throw(ArgumentError("loads must be :design or :frozen, got :$loads"))
    t = cl.fef
    (loads === :frozen || (isempty(t.loads) && isempty(t.sloads) && isempty(t.tloads))) && return cl.F0
    return cl.P - case_Pf(cl, state)
end

"""
    solve(model, state, cl::CaseLoads; solver = nothing, loads = :design) -> Matrix

The pure solve for several load cases at once: one stiffness assembly and
one factorization at the design `state`, one right-hand side per case.
Returns the FULL-space displacements, `n_dof × n_cases` (column `k` is case
`cl.names[k]`; inactive and fixed slots zero).

`loads = :design` (default) re-evaluates every case's element loads at the
state ([`case_Pf`](@ref)): frame line and point loads and self-weight, quad
and triangle shell surface loads. `loads = :frozen` uses the loads at the
model's own state (`cl.F0`). Rigid diaphragms are handled as in the
single-case [`solve`](@ref).

A load combination is a linear combination of the columns: `U * λ` with λ
from [`combination_weights`](@ref). Differentiable with Zygote and
ForwardDiff.
"""
function solve(model::Model{T}, state::ModelState, cl::CaseLoads; solver = nothing,
        loads::Symbol = :design) where {T}
    cache = model.cache::AnalysisCache{T}
    K = assemble_K(cache, state)
    F = _case_F(cl, state, loads)[cache.partition.free, :]
    cm = cache.constraint
    Uf = cm === nothing ? _pure_solve(K, F, solver) : _constrained_solve(cm, K, F, solver)
    return cache.free_embed * Uf
end

"""
    compliance(model, state, cl::CaseLoads; solver = nothing, loads = :design) -> Vector

The external work `Fₖᵀuₖ` of every case at the design `state` (one stiffness
factorization for all of them), `loads` as in [`solve`](@ref)`(model, state, cl)`.

The compliance of a combination `λ` is NOT `λ ⋅ compliance` (it is quadratic
in the loads): it is `(Fλ)ᵀ(Uλ)`, the work of the combined load on the
combined displacement.
"""
function compliance(model::Model{T}, state::ModelState, cl::CaseLoads; solver = nothing,
        loads::Symbol = :design) where {T}
    cache = model.cache::AnalysisCache{T}
    K = assemble_K(cache, state)
    F = _case_F(cl, state, loads)[cache.partition.free, :]
    cm = cache.constraint
    Uf = cm === nothing ? _pure_solve(K, F, solver) : _constrained_solve(cm, K, F, solver)
    return vec(sum(F .* Uf; dims = 1))
end

"""
    combination_weights(cl::CaseLoads, combo::LoadCombination) -> Vector

The factors of a [`LoadCombination`](@ref) as a vector `λ` over the columns
of `cl` (zero for a case the combination does not name), so that `U * λ` and
`F * λ` are the combination's displacements and loads. A case named in the
combination but absent from `cl` throws an `ArgumentError` (a misspelled case
would otherwise drop silently).
"""
function combination_weights(cl::CaseLoads{T}, combo::LoadCombination) where {T}
    for (c, _) in combo.factors
        c in cl.names || throw(ArgumentError("combination :$(combo.name) names case :$c, not in $(cl.names)"))
    end
    f = factor_dict(combo)
    return T[get(f, c, zero(T)) for c in cl.names]
end

"""
    fef_local(cl::CaseLoads, state) -> Matrix   (12 × n_frame_loads)

The LOCAL condensed fixed-end force vectors `q̃` of every frame element load
of the cases at the design `state`, one column per load in `cl.fef` order:
load `j` acts on element `cl.fef.elem[j]` and belongs to case
`cl.names[cl.fcase[j]]`. The multi-case form of
[`fef_local`](@ref)`(cache, state)`: per-case force recovery adds each case's
columns to its own element forces (e.g. slot 7, the end-node axial action, to
the axial force). Pure and differentiable.
"""
function fef_local(cl::CaseLoads, state::ModelState)
    t = cl.fef
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
