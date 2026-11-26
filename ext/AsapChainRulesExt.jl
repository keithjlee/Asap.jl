"""
    AsapChainRulesExt

The custom-rule surface of Asap's differentiable path. Correctness needs
only the first two rules; everything else in the pure pipeline (element
kernels, scatter, embeddings) is mutation-free StaticArrays/SparseArrays
arithmetic that AD engines differentiate natively.

1. `sparse_from_pattern`: constructing a sparse matrix from a frozen
   pattern and a value vector. Reverse: project the (dense or sparse)
   matrix cotangent onto the pattern's stored entries.
2. `solve_free`: the reduced linear solve `u = K⁻¹F`. Reverse: the
   classical adjoint method — one extra back-substitution `λ = K⁻ᵀ ū` on
   the factorization already computed in the forward pass, with
   `K̄ = −λuᵀ` and `F̄ = λ`. K is symmetric, so no re-factorization.

`solve_free` also carries the matching FORWARD rule (`frule`): the
implicit-function theorem `K u̇ = Ḟ − K̇ u`, one extra back-substitution
per tangent on the forward factorization.

The rest are for Zygote: constructor rules for types with custom inner
constructors, a performance rule for the truss kernel, the position seam,
and the linear-cost rules on the per-element maps (end of file; RuleConfig
rules that Enzyme and Mooncake never see).

Loaded automatically as a package extension whenever ChainRulesCore is in
the environment (`[weakdeps]`/`[extensions]` in Project.toml) — zero cost
for users who never differentiate.
"""
module AsapChainRulesExt

using Asap
using ChainRulesCore
using SparseArrays
using LinearAlgebra
using StaticArrays

const HOST = Asap

function ChainRulesCore.rrule(::typeof(HOST.sparse_from_pattern),
    pattern::SparseMatrixCSC, nzval::AbstractVector)

    K = HOST.sparse_from_pattern(pattern, nzval)
    rows = rowvals(pattern)

    function sparse_from_pattern_pullback(ΔK)
        return (NoTangent(), NoTangent(), _pattern_values(unthunk(ΔK), pattern, rows))
    end

    return K, sparse_from_pattern_pullback
end

# extract the pattern-entry cotangents from whichever representation the AD
# engine hands back: a dense/sparse matrix, or a structural Tangent of the
# SparseMatrixCSC (Mooncake's ChainRules bridge) whose nzval matches the
# frozen pattern by construction
function _pattern_values(Δ::AbstractMatrix, pattern::SparseMatrixCSC, rows)
    Δnz = Vector{eltype(pattern)}(undef, nnz(pattern))
    @inbounds for j in 1:size(pattern, 2), p in nzrange(pattern, j)
        Δnz[p] = Δ[rows[p], j]
    end
    return Δnz
end
function _pattern_values(Δ::ChainRulesCore.Tangent, pattern::SparseMatrixCSC, rows)
    nz = Δ.nzval
    nz isa ChainRulesCore.AbstractZero && return zeros(eltype(pattern), nnz(pattern))
    @assert length(nz) == nnz(pattern) "structural sparse cotangent pattern mismatch"
    return collect(nz)
end

# constructor rules: Section/Material have custom inner constructors
# (promotion logic), which Zygote's automatic struct adjoint refuses —
# their pullbacks are simple field pass-throughs
function ChainRulesCore.rrule(::Type{<:Asap.Section}, material::Asap.Material,
    A::Real, Ix::Real, Iy::Real, J::Real)
    sec = Asap.Section(material, A, Ix, Iy, J)
    function Section_pullback(Δ)
        Δs = unthunk(Δ)
        return (NoTangent(), Δs.material, Δs.A, Δs.Ix, Δs.Iy, Δs.J)
    end
    return sec, Section_pullback
end

function ChainRulesCore.rrule(::Type{<:Asap.Material}, E::Real, G::Real, ρ::Real, ν::Real)
    m = Asap.Material(E, G, ρ, ν)
    function Material_pullback(Δ)
        Δm = unthunk(Δ)
        return (NoTangent(), Δm.E, Δm.G, Δm.ρ, Δm.ν)
    end
    return m, Material_pullback
end

# constructor rules for the end-condition types (custom inner constructors
# with promotion/assertions block Zygote's automatic struct adjoint) —
# needed when connection stiffness is a design variable
function ChainRulesCore.rrule(::Type{<:Asap.EndSprings}, kx::Real, kt::Real, ky::Real, kz::Real)
    e = Asap.EndSprings(kx, kt, ky, kz)
    function EndSprings_pullback(Δ)
        Δe = unthunk(Δ)
        return (NoTangent(), Δe.kx, Δe.kt, Δe.ky, Δe.kz)
    end
    return e, EndSprings_pullback
end

function ChainRulesCore.rrule(::Type{<:Asap.EndConditions}, e1::Asap.EndSprings, e2::Asap.EndSprings)
    ec = Asap.EndConditions(e1, e2)
    function EndConditions_pullback(Δ)
        Δc = unthunk(Δ)
        return (NoTangent(), Δc.e1, Δc.e2)
    end
    return ec, EndConditions_pullback
end

# PERFORMANCE rule (correctness needs nothing here — Zygote traverses the
# kernel natively, verified by tests — but the analytic pullback is ~50×
# cheaper and truss kernels dominate large sizing problems):
#   B(section, x1, x2) = (EA/L)·[nnᵀ −nnᵀ; −nnᵀ nnᵀ],  v = x2−x1, n = v/L
# Given block-combined symmetrized cotangent G = sym(ΔB₁₁ − ΔB₁₂ − ΔB₂₁ + ΔB₂₂):
#   ΔEA = (nᵀGn)/L,   Δv = (2c/L)·G·n − (3c/L)·(nᵀGn)·n   with c = EA/L
function ChainRulesCore.rrule(::typeof(Asap.truss_stiffness),
    section::Asap.AbstractSection, x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real})

    v = SVector{3}(x2) - SVector{3}(x1)
    L = norm(v)
    n = v / L
    ea = Asap.EA(section)
    c = ea / L
    B = Asap.truss_stiffness(section, x1, x2)

    function truss_stiffness_pullback(ΔB)
        Δ = unthunk(ΔB)
        G11 = SMatrix{3,3}(ntuple(k -> Δ[(k-1)%3+1, (k-1)÷3+1], 9))
        G12 = SMatrix{3,3}(ntuple(k -> Δ[(k-1)%3+1, (k-1)÷3+4], 9))
        G21 = SMatrix{3,3}(ntuple(k -> Δ[(k-1)%3+4, (k-1)÷3+1], 9))
        G22 = SMatrix{3,3}(ntuple(k -> Δ[(k-1)%3+4, (k-1)÷3+4], 9))
        G = G11 - G12 - G21 + G22
        Gs = (G + G') / 2

        nGn = dot(n, Gs, n)
        ΔEA = nGn / L
        Δv = (2c / L) * (Gs * n) - (3c / L) * nGn * n

        return (NoTangent(), _section_ea_tangent(section, ΔEA), -Δv, Δv)
    end

    return B, truss_stiffness_pullback
end

_section_ea_tangent(s::Asap.Section, ΔEA) =
    ChainRulesCore.Tangent{typeof(s)}(
        material=ChainRulesCore.Tangent{typeof(s.material)}(E=ΔEA * s.A),
        A=ΔEA * s.material.E)
_section_ea_tangent(s::Asap.RigiditySection, ΔEA) =
    ChainRulesCore.Tangent{typeof(s)}(EA=ΔEA)

# boundary rule between the plain position matrix (AD input container) and
# the static vectors the kernels use internally — avoids the known
# Zygote/StaticArrays tangent-projection mismatch at exactly this seam
function ChainRulesCore.rrule(::typeof(HOST._position), X::AbstractMatrix, i::Int)
    y = HOST._position(X, i)
    function _position_pullback(Δ)
        Δv = _flatten3(unthunk(Δ))
        ΔX = ChainRulesCore.@thunk begin
            dX = zeros(eltype(X), size(X))
            dX[1, i] = Δv[1]
            dX[2, i] = Δv[2]
            dX[3, i] = Δv[3]
            dX
        end
        return (NoTangent(), ΔX, NoTangent())
    end
    return y, _position_pullback
end

# static-vector cotangents arrive in several structural wrappings depending
# on the AD engine's path (natural SVector, Tangent{SVector}(data = ...),
# nested Tangent{Any, Tuple}); flatten them all to a plain 3-tuple
_flatten3(Δ::AbstractVector) = (Δ[1], Δ[2], Δ[3])
_flatten3(Δ::NTuple{3,Any}) = Δ
_flatten3(Δ::ChainRulesCore.Tangent) = _flatten3(ChainRulesCore.backing(Δ))
_flatten3(Δ::NamedTuple) = _flatten3(first(values(Δ)))

# The K-cotangent is returned PROJECTED onto K's sparsity pattern (the
# nonzero entries of −λuᵀ that K actually stores). In Asap's pipeline K
# always comes from sparse_from_pattern, whose pullback reads exactly
# those entries — so this loses nothing, avoids materializing the dense
# n×n outer product, and keeps the cotangent structurally compatible with
# sparse-tangent AD engines (Mooncake).
function _project_K_cotangent(K::SparseMatrixCSC, λ, u)
    rows = rowvals(K)
    vals = Vector{eltype(K)}(undef, nnz(K))
    @inbounds for j in 1:size(K, 2), p in nzrange(K, j)
        vals[p] = -_rowdot(λ, u, rows[p], j)
    end
    return SparseMatrixCSC(size(K, 1), size(K, 2), copy(K.colptr), copy(K.rowval), vals)
end

# NB two Enzyme-@import_rrule constraints on the 2-arg rule (the imported
# one): the pullback must be FLAT (literal return tuple — no shared inner
# closure with a splatted return), and the factorization must be the RAW
# concretely-typed object from the one-arg _factorize — the
# FactorizationCache wrapper's Any-typed fields put dynamic dispatch
# inside the imported rule ("Illegal calling convention fixup").
function ChainRulesCore.rrule(::typeof(HOST.solve_free),
    K::SparseMatrixCSC, F::AbstractVecOrMat)
    fact = HOST._factorize(K)
    u = fact \ F
    function solve_free_pullback(Δu)
        λ = fact \ collect(unthunk(Δu))  # K symmetric: Kᵀ = K, reuse the factorization
        return (NoTangent(), @thunk(_project_K_cotangent(K, λ, u)), λ)
    end
    return u, solve_free_pullback
end

function ChainRulesCore.rrule(::typeof(HOST.solve_free), solver,
    K::SparseMatrixCSC, F::AbstractVecOrMat)
    fact = HOST._factorize(solver, K)
    u = fact \ F
    function solve_free_solver_pullback(Δu)
        λ = fact \ collect(unthunk(Δu))
        return (NoTangent(), NoTangent(), @thunk(_project_K_cotangent(K, λ, u)), λ)
    end
    return u, solve_free_solver_pullback
end

@inline _rowdot(λ::AbstractVector, u::AbstractVector, i::Int, j::Int) = λ[i] * u[j]
@inline _rowdot(λ::AbstractMatrix, u::AbstractMatrix, i::Int, j::Int) =
    dot(view(λ, i, :), view(u, j, :))

"""
Forward rule for the reduced solve — the implicit-function theorem:
`K u̇ = Ḟ − K̇ u`, one extra back-substitution on the forward
factorization per tangent. Consumed by forward-mode ChainRules engines
and bridged to Enzyme's forward mode via `Enzyme.@import_frule`
(`AsapEnzymeExt`). ForwardDiff does NOT read frules — its Dual path is
`AsapForwardDiffExt`.
"""
function _solve_free_forward(solver, K, F, ΔK, ΔF)
    fact = HOST._factorize(solver, K)
    u = fact \ F
    ΔKm = unthunk(ΔK)
    ΔFm = unthunk(ΔF)
    rhs = ΔFm isa AbstractZero ? zero(u) : collect(ΔFm)
    ΔKm isa AbstractZero || (rhs = rhs .- ΔKm * u)
    return u, fact \ rhs
end

ChainRulesCore.frule((_, ΔK, ΔF), ::typeof(HOST.solve_free),
    K::SparseMatrixCSC, F::AbstractVecOrMat) =
    _solve_free_forward(nothing, K, F, ΔK, ΔF)

ChainRulesCore.frule((_, _, ΔK, ΔF), ::typeof(HOST.solve_free), solver,
    K::SparseMatrixCSC, F::AbstractVecOrMat) =
    _solve_free_forward(solver, K, F, ΔK, ΔF)

# ─────────────────────────────────────────────────────────────────────────────
# Linear-cost reverse rules for the per-element maps (Zygote; any RuleConfig
# reverse engine)
# ─────────────────────────────────────────────────────────────────────────────
# The pure path evaluates one closure per element that indexes the captured
# state arrays (`_position(state.X, i)`, `state.sections[e]`, `state.ends[e]`,
# `t.loads[j]`). Under Zygote each such index pulls back into a FULL-size
# cotangent array, so one gradient costs O(n²) in the number of elements.
#
# These rules leave the primal code exactly as it is and, in the pullback,
# differentiate each element's kernel on its own gathered arguments
# (`rrule_via_ad`, so the kernel rules of AsapFrameADExt/AsapShellADExt apply),
# accumulating straight into ONE cotangent per state field.
#
# They take a RuleConfig on purpose: Enzyme (`@import_rrule`) and Mooncake
# (`@from_rrule`) only ever import config-free rules for named kernels, so they
# keep differentiating the unchanged primal natively. Restructuring the primal
# itself (gather + `map` over `Vector{Any}`) fixed Zygote but made Enzyme
# reverse throw and Enzyme forward silently wrong.

const _RevCfg = RuleConfig{>:HasReverseMode}

# Cotangent of `state` from per-element units. `item(j)` returns
# `(f, s, nodes, e, ec, m)`: the kernel `f(section, positions::Tuple, ends)`
# (capturing only small constants), the section index `s` into
# `state.sections`, the node indices, the index `e` into `state.ends` (0 when
# the ends are the model's constants), the ends value, and the block length.
function _itemwise_tangent(config, state::ModelState, n::Int, item, Δ::AbstractVector)
    dX = zero(state.X)
    dS = Vector{Any}(undef, length(state.sections))
    fill!(dS, ZeroTangent())
    dE = state.ends === nothing ? ZeroTangent() : fill!(Vector{Any}(undef, length(state.ends)), ZeroTangent())
    off = 0
    for j in 1:n
        f, s, nodes, e, ec, m = item(j)
        xs = map(i -> HOST._position(state.X, i), nodes)
        _, pb = rrule_via_ad(config, f, state.sections[s], xs, ec)
        _, Δs, Δx, Δe = pb(Δ[off+1:off+m])
        off += m
        dS[s] = _add(dS[s], Δs)
        _accum_positions!(dX, nodes, Δx)
        e > 0 && (dE[e] = _add(dE[e], Δe))
    end
    return Tangent{typeof(state)}(X = dX, sections = dS, ends = dE)
end

_add(a, b) = a + _natural(unthunk(b))
_add(a, ::AbstractZero) = a

# `rrule_via_ad` hands static-array fields of a struct cotangent back as
# structural `Tangent{SMatrix}(data = …)`; user code upstream (e.g. `S .* fS`
# in a section constructor) needs the natural array, as plain Zygote gives it
_natural(Δ) = Δ
_natural(Δ::Tangent{P}) where {P<:StaticArray} =
    (d = ChainRulesCore.backing(Δ).data; d isa AbstractZero ? d : P(map(unthunk, d)))
function _natural(Δ::Tangent{P,<:NamedTuple}) where {P}
    nt = map(x -> _natural(unthunk(x)), ChainRulesCore.backing(Δ))
    return Tangent{P,typeof(nt)}(nt)
end

_accum_positions!(dX, nodes, ::AbstractZero) = dX
_accum_positions!(dX, nodes, Δx::Tangent) = _accum_positions!(dX, nodes, ChainRulesCore.backing(Δx))
function _accum_positions!(dX, nodes, Δx::Tuple)
    for (k, i) in enumerate(nodes)
        Δk = unthunk(Δx[k])
        Δk isa AbstractZero && continue
        v = _flatten3(Δk)
        dX[1, i] += v[1]; dX[2, i] += v[2]; dX[3, i] += v[3]
    end
    return dX
end

_cotvec(Δ) = (d = unthunk(Δ); d isa AbstractZero ? d : collect(vec(d)))

const _F64State = ModelState{<:Base.IEEEFloat}

# one frame's active stiffness entries from plain arguments: the same two lines
# as `Asap._group_entries` for frames. Deliberately NOT factored out of the
# primal: Enzyme's type analysis failed on that helper (frame_joints reverse,
# "Cannot deduce type of load"), so the primal stays exactly as Enzyme knows it
function _frame_entries(section, ec, x1, x2, rollangle, slots)
    Ke = Asap.frame_stiffness(section, ec, x1, x2, rollangle)
    return length(slots) == 12 ? vec(Ke) : vec(Ke[slots, slots])
end

# frame stiffness entries of one group
function ChainRulesCore.rrule(config::_RevCfg, ::typeof(HOST._group_entries_all),
        g::HOST.ElementGroup{<:Asap.FrameElement}, state::_F64State)
    V = HOST._group_entries_all(g, state)
    own = state.ends === nothing
    function item(e)
        mi, θ, slots = g.model_indices[e], g.rollangles[e], g.slots[e]
        f = (sec, xs, ec) -> _frame_entries(sec, ec, xs[1], xs[2], θ, slots)
        return f, mi, (g.i1[e], g.i2[e]), own ? 0 : mi, own ? g.endss[e] : state.ends[mi], length(slots)^2
    end
    function frame_group_pullback(ΔV)
        Δ = _cotvec(ΔV)
        Δ isa AbstractZero && return (NoTangent(), NoTangent(), ZeroTangent())
        return (NoTangent(), NoTangent(), _itemwise_tangent(config, state, length(g.elements), item, Δ))
    end
    return V, frame_group_pullback
end

# shell stiffness entries of one group (quadrilaterals N = 4, triangles N = 3)
function ChainRulesCore.rrule(config::_RevCfg, ::typeof(HOST._group_entries_nnode_all),
        g::HOST.ElementGroup{<:Asap.ShellElement}, state::_F64State, v::Val{N}) where {N}
    V = HOST._group_entries_nnode_all(g, state, v)
    function item(e)
        q = HOST._group_quadedge(g, e, Val(N))
        shellK = (sec, xs, ec) -> vec(HOST._shell_stiffness_q(sec, xs, q))
        return (shellK, g.model_indices[e], ntuple(k -> g.inodes[k, e], Val(N)), 0, nothing, (6N)^2)
    end
    function shell_group_pullback(ΔV)
        Δ = _cotvec(ΔV)
        Δ isa AbstractZero && return (NoTangent(), NoTangent(), ZeroTangent(), NoTangent())
        return (NoTangent(), NoTangent(), _itemwise_tangent(config, state, length(g.elements), item, Δ), NoTangent())
    end
    return V, shell_group_pullback
end

# design-dependent frame element loads (global, scattered) and their local q̃
function _frame_fef_item(t, state)
    own = state.ends === nothing
    return function (j)
        load, e, θ = t.loads[j], t.elem[j], t.rollangle[j]
        return load, e, (t.i1[j], t.i2[j]), own ? 0 : e, own ? t.ends[j] : state.ends[e], θ
    end
end

function ChainRulesCore.rrule(config::_RevCfg, ::typeof(HOST._frame_fef),
        t::HOST.FEFTable, state::_F64State)
    y = HOST._frame_fef(t, state)
    base = _frame_fef_item(t, state)
    function item(j)
        load, s, nodes, e, ec, θ = base(j)
        f = (sec, xs, ec) -> Asap.element_fef_global(load, sec, ec, xs[1], xs[2], θ)
        return f, s, nodes, e, ec, 12
    end
    function frame_fef_pullback(Δy)
        Δ = _cotvec(Δy)
        Δ isa AbstractZero && return (NoTangent(), NoTangent(), ZeroTangent())
        return (NoTangent(), NoTangent(), _itemwise_tangent(config, state, length(t.loads), item, t.scatter' * Δ))
    end
    return y, frame_fef_pullback
end

function ChainRulesCore.rrule(config::_RevCfg, ::typeof(Asap.fef_local),
        cache::HOST.AnalysisCache, state::_F64State)
    y = Asap.fef_local(cache, state)
    t = cache.fef
    base = _frame_fef_item(t, state)
    function item(j)
        load, s, nodes, e, ec, θ = base(j)
        f = (sec, xs, ec) -> Asap.condensed_fef(load, sec, ec, xs[1], xs[2], θ)
        return f, s, nodes, e, ec, 12
    end
    function fef_local_pullback(Δy)
        Δ = _cotvec(Δy)
        (Δ isa AbstractZero || isempty(t.loads)) && return (NoTangent(), NoTangent(), ZeroTangent())
        return (NoTangent(), NoTangent(), _itemwise_tangent(config, state, length(t.loads), item, Δ))
    end
    return y, fef_local_pullback
end

# design-dependent shell surface loads
function ChainRulesCore.rrule(config::_RevCfg, ::typeof(HOST._shell_fef),
        t::HOST.FEFTable, state::_F64State)
    y = HOST._shell_fef(t, state)
    function item(j)
        load = t.sloads[j]
        f = (sec, xs, ec) -> HOST._shell_fef_global(load, sec, xs)
        return f, t.selem[j], ntuple(k -> t.snodes[k, j], 4), 0, nothing, 24
    end
    function shell_fef_pullback(Δy)
        Δ = _cotvec(Δy)
        Δ isa AbstractZero && return (NoTangent(), NoTangent(), ZeroTangent())
        return (NoTangent(), NoTangent(), _itemwise_tangent(config, state, length(t.sloads), item, t.sscatter' * Δ))
    end
    return y, shell_fef_pullback
end

# design-dependent triangular shell surface loads
function ChainRulesCore.rrule(config::_RevCfg, ::typeof(HOST._shell3_fef),
        t::HOST.FEFTable, state::_F64State)
    y = HOST._shell3_fef(t, state)
    function item(j)
        load = t.tloads[j]
        f = (sec, xs, ec) -> HOST._shell_fef_global(load, sec, xs)
        return f, t.telem[j], ntuple(k -> t.tnodes[k, j], Val(3)), 0, nothing, 18
    end
    function shell3_fef_pullback(Δy)
        Δ = _cotvec(Δy)
        Δ isa AbstractZero && return (NoTangent(), NoTangent(), ZeroTangent())
        return (NoTangent(), NoTangent(), _itemwise_tangent(config, state, length(t.tloads), item, t.tscatter' * Δ))
    end
    return y, shell3_fef_pullback
end

# the stacked (unscattered) fixed-end vectors of the load cases (`CaseLoads`):
# the same per-load kernels as above, the pullback takes the stacked cotangent
# directly (the case scatter is a plain sparse matvec Zygote handles itself)
function ChainRulesCore.rrule(config::_RevCfg, ::typeof(HOST._frame_fef_q),
        t::HOST.FEFTable, state::_F64State)
    y = HOST._frame_fef_q(t, state)
    base = _frame_fef_item(t, state)
    function item(j)
        load, s, nodes, e, ec, θ = base(j)
        f = (sec, xs, ec) -> Asap.element_fef_global(load, sec, ec, xs[1], xs[2], θ)
        return f, s, nodes, e, ec, 12
    end
    function frame_fef_q_pullback(Δy)
        Δ = _cotvec(Δy)
        Δ isa AbstractZero && return (NoTangent(), NoTangent(), ZeroTangent())
        return (NoTangent(), NoTangent(), _itemwise_tangent(config, state, length(t.loads), item, Δ))
    end
    return y, frame_fef_q_pullback
end

function ChainRulesCore.rrule(config::_RevCfg, ::typeof(HOST._shell_fef_q),
        t::HOST.FEFTable, state::_F64State)
    y = HOST._shell_fef_q(t, state)
    function item(j)
        load = t.sloads[j]
        f = (sec, xs, ec) -> HOST._shell_fef_global(load, sec, xs)
        return f, t.selem[j], ntuple(k -> t.snodes[k, j], 4), 0, nothing, 24
    end
    function shell_fef_q_pullback(Δy)
        Δ = _cotvec(Δy)
        Δ isa AbstractZero && return (NoTangent(), NoTangent(), ZeroTangent())
        return (NoTangent(), NoTangent(), _itemwise_tangent(config, state, length(t.sloads), item, Δ))
    end
    return y, shell_fef_q_pullback
end

function ChainRulesCore.rrule(config::_RevCfg, ::typeof(HOST._shell3_fef_q),
        t::HOST.FEFTable, state::_F64State)
    y = HOST._shell3_fef_q(t, state)
    function item(j)
        load = t.tloads[j]
        f = (sec, xs, ec) -> HOST._shell_fef_global(load, sec, xs)
        return f, t.telem[j], ntuple(k -> t.tnodes[k, j], Val(3)), 0, nothing, 18
    end
    function shell3_fef_q_pullback(Δy)
        Δ = _cotvec(Δy)
        Δ isa AbstractZero && return (NoTangent(), NoTangent(), ZeroTangent())
        return (NoTangent(), NoTangent(), _itemwise_tangent(config, state, length(t.tloads), item, Δ))
    end
    return y, shell3_fef_q_pullback
end

function ChainRulesCore.rrule(config::_RevCfg, ::typeof(Asap.fef_local),
        cl::Asap.CaseLoads, state::_F64State)
    y = Asap.fef_local(cl, state)
    t = cl.fef
    base = _frame_fef_item(t, state)
    function item(j)
        load, s, nodes, e, ec, θ = base(j)
        f = (sec, xs, ec) -> Asap.condensed_fef(load, sec, ec, xs[1], xs[2], θ)
        return f, s, nodes, e, ec, 12
    end
    function case_fef_local_pullback(Δy)
        Δ = _cotvec(Δy)
        (Δ isa AbstractZero || isempty(t.loads)) && return (NoTangent(), NoTangent(), ZeroTangent())
        return (NoTangent(), NoTangent(), _itemwise_tangent(config, state, length(t.loads), item, Δ))
    end
    return y, case_fef_local_pullback
end

end # module
