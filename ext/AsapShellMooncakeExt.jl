"""
    AsapShellMooncakeExt

Mooncake bridge for the shell kernels: Mooncake's own rule derivation
overflows the stack on the 24×24 static-matrix shell kernels, so the
ChainRules rules of `AsapShellADExt` (local ForwardDiff Jacobians) are
imported with `Mooncake.@from_rrule` — the same mechanism `AsapMooncakeExt`
uses for the sparse construction and the linear solve.

Active when Mooncake, ChainRulesCore and ForwardDiff are all loaded.
"""
module AsapShellMooncakeExt

using Asap, Mooncake, ChainRulesCore, ForwardDiff, StaticArrays

const M3 = SMatrix{3,3,Float64,9}
const M2 = SMatrix{2,2,Float64,4}
const V3 = SVector{3,Float64}

# Bridge adapter: the rules return static-array cotangents; Mooncake stores a
# static array's rdata as `(data = NTuple,)` and lacks the ChainRules
# conversion for it — accumulate the tuple
function Mooncake.increment_and_get_rdata!(::Mooncake.NoFData,
    r::Mooncake.RData{NamedTuple{(:data,),Tuple{NTuple{N,Float64}}}},
    t::StaticArray{S,Float64}) where {N,S}
    return Mooncake.RData((data = map(+, r.data.data, Tuple(t)),))
end

Mooncake.@from_rrule(
    Mooncake.DefaultCtx,
    Tuple{typeof(Asap._shell_K_plain),M3,M3,M2,Float64,V3,V3,V3,V3},
)

Mooncake.@from_rrule(
    Mooncake.DefaultCtx,
    Tuple{typeof(Asap._shell_fef_plain),M3,V3,V3,V3,V3,V3},
)

end # module
