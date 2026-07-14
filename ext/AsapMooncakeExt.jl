"""
    AsapMooncakeExt

Mooncake bridge: Mooncake does not consume ChainRulesCore rules
automatically, so the handful of rules that make Asap's pure path
differentiable (the sparse construction and the CHOLMOD-backed linear
solve — foreign calls no AD engine can traverse) are imported explicitly
with `Mooncake.@from_rrule`. Everything else in the pipeline Mooncake
differentiates natively.

Activates automatically when Mooncake is loaded alongside Asap.
"""
module AsapMooncakeExt

using Asap
using Mooncake
using SparseArrays

# A CHOLMOD factorization is opaque foreign solver state — never
# differentiable. Without this declaration Mooncake tries to build tangent
# storage for a Factor captured in a solved model's cache (via the objective
# closure) and fails deep in its internals.
Mooncake.tangent_type(::Type{<:SparseArrays.CHOLMOD.Factor}) = Mooncake.NoTangent

# Bridge adapter: Asap's solve_free rule returns its K-cotangent as a
# SparseMatrixCSC sharing the primal's pattern; Mooncake's ChainRules
# interop lacks a conversion for sparse cotangents onto its structural
# sparse fdata — supply it (accumulate stored values; patterns match by
# construction in Asap's rules).
function Mooncake.increment_and_get_rdata!(
    f::Mooncake.FData{<:NamedTuple{(:m, :n, :colptr, :rowval, :nzval)}},
    ::Mooncake.NoRData,
    t::SparseMatrixCSC{Float64,Int},
)
    @assert length(f.data.nzval) == SparseArrays.nnz(t) "sparse cotangent pattern mismatch"
    f.data.nzval .+= t.nzval
    return Mooncake.NoRData()
end

Mooncake.@from_rrule(
    Mooncake.DefaultCtx,
    Tuple{typeof(Asap.solve_free),SparseMatrixCSC{Float64,Int},Vector{Float64}},
)

Mooncake.@from_rrule(
    Mooncake.DefaultCtx,
    Tuple{typeof(Asap.solve_free),SparseMatrixCSC{Float64,Int},Matrix{Float64}},
)

Mooncake.@from_rrule(
    Mooncake.DefaultCtx,
    Tuple{typeof(Asap.sparse_from_pattern),SparseMatrixCSC{Float64,Int},Vector{Float64}},
)

end # module
