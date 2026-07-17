"""
    AsapForwardDiffExt

ForwardDiff bridge: `Dual`-specialized methods for the one operation Dual
numbers cannot pass through — the CHOLMOD/LinearSolve-backed sparse solve
(`solve_free`), whose factorization is Float64-only foreign code.

The methods implement the implicit-function theorem directly:

    K(θ) u(θ) = F(θ)   ⇒   K u̇ᵏ = Ḟᵏ − K̇ᵏ u

so the VALUE system is factorized once, and every partial is one extra
multi-RHS back-substitution on that same factorization. This is what makes
forward-mode gradients scale with (design variables ÷ chunk size) cheap
extra solves — the economics that favor forward mode whenever design
variables are fewer than constraint outputs.

Everything else in Asap's pure pipeline (assembly kernels, section
promotion, geometry) is generic Julia code that Dual numbers traverse
natively — no further rules needed.

Activates automatically when ForwardDiff is loaded alongside Asap:

    using DifferentiationInterface, ForwardDiff
    DifferentiationInterface.gradient(obj, AutoForwardDiff(), x0)
"""
module AsapForwardDiffExt

using Asap
using ForwardDiff
using ForwardDiff: Dual, Partials, value, partials
using SparseArrays

# value/partial extraction that tolerates plain (non-Dual) arrays: in mixed
# systems (Dual K, constant F — the standard case, loads are design-constant)
# the missing partials are zero, expressed as the broadcast-neutral `false`
_values(A::AbstractArray{<:Dual}) = value.(A)
_values(A::AbstractArray) = A
_partials(A::AbstractArray{<:Dual}, k::Int) = partials.(A, k)
_partials(::AbstractArray, ::Int) = false

_same_pattern(K::SparseMatrixCSC, nz::AbstractVector) =
    SparseMatrixCSC(size(K, 1), size(K, 2), K.colptr, K.rowval, nz)

function _solve_free_dual(::Type{Dual{Tg,V,N}}, K::SparseMatrixCSC,
    F::AbstractVecOrMat) where {Tg,V,N}

    Kv = _same_pattern(K, _values(K.nzval))
    fact = Asap._factorize(Kv)               # one factorization for value + all N partials
    u = fact \ _values(F)
    K_dual = eltype(K) <: Dual
    dus = ntuple(N) do k
        rhs = K_dual ?
              _partials(F, k) .- _same_pattern(K, _partials(K.nzval, k)) * u :
              _partials(F, k)                # plain K ⇒ F carries the Duals
        fact \ rhs
    end
    return Dual{Tg}.(u, Partials.(tuple.(dus...)))
end

# Dual stiffness (positions/sections/joints carry derivatives). Concrete
# vector/matrix methods mirror the base signatures — Dual <: Real, so a
# single loose method would be ambiguous against Asap's own.
Asap.solve_free(K::SparseMatrixCSC{D}, F::AbstractVector) where {D<:Dual} =
    _solve_free_dual(D, K, F)
Asap.solve_free(K::SparseMatrixCSC{D}, F::AbstractMatrix) where {D<:Dual} =
    _solve_free_dual(D, K, F)

# plain stiffness, Dual RHS (load-only derivatives): u̇ᵏ = K⁻¹Ḟᵏ
Asap.solve_free(K::SparseMatrixCSC{<:AbstractFloat}, F::AbstractVector{D}) where {D<:Dual} =
    _solve_free_dual(D, K, F)
Asap.solve_free(K::SparseMatrixCSC{<:AbstractFloat}, F::AbstractMatrix{D}) where {D<:Dual} =
    _solve_free_dual(D, K, F)

end # module
