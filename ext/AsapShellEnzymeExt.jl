"""
    AsapShellEnzymeExt

Enzyme bridge for the shell kernels. Differentiating them natively fails in
Enzyme's type analysis (the 24×24 static-matrix code exceeds its register
size limits), so the ChainRules rules of `AsapShellADExt` (local ForwardDiff
Jacobians) are imported with `Enzyme.@import_rrule` — the same mechanism
`AsapEnzymeExt` uses for the solve. REVERSE mode only (no forward rules are
imported for the shell kernels; use ForwardDiff for forward mode).

Active when Enzyme, ChainRulesCore and ForwardDiff are all loaded. Use the
usual Asap Enzyme settings: `set_runtime_activity(Reverse)` and a `Const`
function annotation.
"""
module AsapShellEnzymeExt

using Asap, Enzyme, ChainRulesCore, ForwardDiff, StaticArrays

# registered at load time: Enzyme's @import_rrule lives in ITS ChainRulesCore
# extension, which is not reliably active while this extension precompiles
function __init__()
    @eval begin
        Enzyme.@import_rrule(typeof(Asap._shell_K_plain), SMatrix{3,3,Float64,9}, SMatrix{3,3,Float64,9},
            SMatrix{2,2,Float64,4}, Float64, SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64})
        Enzyme.@import_rrule(typeof(Asap._shell_fef_plain), SMatrix{3,3,Float64,9}, SVector{3,Float64},
            SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64})
    end
    return nothing
end

end # module
