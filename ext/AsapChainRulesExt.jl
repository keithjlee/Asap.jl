"""
    AsapChainRulesExt

The complete custom-rule surface of Asap's differentiable path — exactly
two reverse rules. Everything else in the pure pipeline (element kernels,
scatter, embeddings) is mutation-free StaticArrays/SparseArrays arithmetic
that AD engines differentiate natively.

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

function ChainRulesCore.rrule(::typeof(HOST.solve_free),
    K::SparseMatrixCSC, F::AbstractVecOrMat)

    fact = HOST._factorize(K)
    u = fact \ F

    # The K-cotangent is returned PROJECTED onto K's sparsity pattern (the
    # nonzero entries of −λuᵀ that K actually stores). In Asap's pipeline K
    # always comes from sparse_from_pattern, whose pullback reads exactly
    # those entries — so this loses nothing, avoids materializing the dense
    # n×n outer product, and keeps the cotangent structurally compatible
    # with sparse-tangent AD engines (Mooncake).
    function solve_free_pullback(Δu)
        λ = fact \ collect(unthunk(Δu))  # K symmetric: Kᵀ = K, reuse the factorization
        ΔK = @thunk begin
            rows = rowvals(K)
            vals = Vector{eltype(K)}(undef, nnz(K))
            @inbounds for j in 1:size(K, 2), p in nzrange(K, j)
                vals[p] = -_rowdot(λ, u, rows[p], j)
            end
            SparseMatrixCSC(size(K, 1), size(K, 2), copy(K.colptr), copy(K.rowval), vals)
        end
        return (NoTangent(), ΔK, λ)
    end

    return u, solve_free_pullback
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
function ChainRulesCore.frule((_, ΔK, ΔF), ::typeof(HOST.solve_free),
    K::SparseMatrixCSC, F::AbstractVecOrMat)

    fact = HOST._factorize(K)
    u = fact \ F
    ΔKm = unthunk(ΔK)
    ΔFm = unthunk(ΔF)
    rhs = ΔFm isa AbstractZero ? zero(u) : collect(ΔFm)
    ΔKm isa AbstractZero || (rhs = rhs .- ΔKm * u)
    return u, fact \ rhs
end

end # module
