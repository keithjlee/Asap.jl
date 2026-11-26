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
    # never while precompiling a DEPENDENT package (e.g. AsapShellEnzymeExt loads this
    # one during its own precompilation): evaluating into a closed module then
    # breaks incremental compilation (a hard error on Julia 1.12). The rules are
    # registered when the extension is actually loaded.
    ccall(:jl_generating_output, Cint, ()) == 1 && return nothing
    @eval begin
        Enzyme.@import_rrule(typeof(Asap._shell_K_plain), SMatrix{3,3,Float64,9}, SMatrix{3,3,Float64,9},
            SMatrix{2,2,Float64,4}, Float64, SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64})
        Enzyme.@import_rrule(typeof(Asap._shell_fef_plain), SMatrix{3,3,Float64,9}, SVector{3,Float64},
            SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64})
        # triangles (3 positions)
        Enzyme.@import_rrule(typeof(Asap._shell_K_plain), SMatrix{3,3,Float64,9}, SMatrix{3,3,Float64,9},
            SMatrix{2,2,Float64,4}, Float64, SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64})
        Enzyme.@import_rrule(typeof(Asap._shell_fef_plain), SMatrix{3,3,Float64,9}, SVector{3,Float64},
            SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64})
        # triangles beside a quadrilateral (per-edge drilling lumping)
        Enzyme.@import_rrule(typeof(Asap._shell3_K_plain), SMatrix{3,3,Float64,9}, SMatrix{3,3,Float64,9},
            SMatrix{2,2,Float64,4}, Float64, SVector{3,Float64}, SVector{3,Float64}, SVector{3,Float64},
            NTuple{3,Bool})
    end
    return nothing
end

end # module
