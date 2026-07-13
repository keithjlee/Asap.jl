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

Loaded automatically as a package extension whenever ChainRulesCore is in
the environment (`[weakdeps]`/`[extensions]` in Project.toml) — zero cost
for users who never differentiate.
"""
module AsapChainRulesExt

using Asap
using ChainRulesCore
using SparseArrays
using LinearAlgebra

const HOST = Asap

function ChainRulesCore.rrule(::typeof(HOST.sparse_from_pattern),
    pattern::SparseMatrixCSC, nzval::AbstractVector)

    K = HOST.sparse_from_pattern(pattern, nzval)
    rows = rowvals(pattern)

    function sparse_from_pattern_pullback(ΔK)
        ΔKd = unthunk(ΔK)
        Δnz = similar(nzval)
        @inbounds for j in 1:size(pattern, 2), p in nzrange(pattern, j)
            Δnz[p] = ΔKd[rows[p], j]
        end
        return (NoTangent(), NoTangent(), Δnz)
    end

    return K, sparse_from_pattern_pullback
end

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
    K::SparseMatrixCSC, F::AbstractVector)

    fact = HOST._factorize(K)
    u = fact \ F

    function solve_free_pullback(Δu)
        λ = fact \ unthunk(Δu)          # K symmetric: Kᵀ = K, reuse the factorization
        ΔK = @thunk(-λ * u')
        return (NoTangent(), ΔK, λ)
    end

    return u, solve_free_pullback
end

end # module
