"""
    AsapEnzymeExt

Enzyme bridge: imports Asap's ChainRules rules for the two operations no
LLVM-level AD can traverse (the frozen-pattern sparse construction and the
CHOLMOD-backed linear solve) via `Enzyme.@import_rrule`. Everything else in
the pure pipeline Enzyme differentiates natively.

USAGE NOTE: Asap's differentiable path captures the (constant) model inside
closures. Enzyme needs two annotations for this pattern:
runtime activity on the mode, and `Const` treatment of the objective
closure (its captured params are data, not differentiable state):

    using DifferentiationInterface, Enzyme
    backend = AutoEnzyme(;
        mode = Enzyme.set_runtime_activity(Enzyme.Reverse),
        function_annotation = Enzyme.Const)
    DifferentiationInterface.gradient(x -> compliance(solve_structure(x, p), p), backend, x0)

Forward mode works the same way with `Enzyme.Forward` (the solve_free
frule is imported alongside the rrules).

LOADING NOTE: this extension requires BOTH Enzyme and ChainRulesCore in
the session (`[extensions]` condition). Zygote/ChainRules-based setups
load ChainRulesCore implicitly; with Enzyme alone add
`using ChainRulesCore` — otherwise no rules register and Enzyme crashes
inside CHOLMOD with an IllegalTypeAnalysisException.
"""
module AsapEnzymeExt

using Asap
using Enzyme
using SparseArrays

# Enzyme's @import_rrule is implemented inside ITS ChainRulesCore extension,
# which is not reliably active while THIS extension precompiles — register
# the imported rules at load time instead.
function __init__()
    @eval begin
        Enzyme.@import_rrule(typeof(Asap.solve_free), SparseMatrixCSC{Float64,Int}, Vector{Float64})
        Enzyme.@import_rrule(typeof(Asap.solve_free), SparseMatrixCSC{Float64,Int}, Matrix{Float64})
        Enzyme.@import_rrule(typeof(Asap.sparse_from_pattern), SparseMatrixCSC{Float64,Int}, Vector{Float64})
        # forward mode: the solve_free frule (implicit-function theorem) —
        # everything else Enzyme's forward pass traverses natively
        Enzyme.@import_frule(typeof(Asap.solve_free), SparseMatrixCSC{Float64,Int}, Vector{Float64})
        Enzyme.@import_frule(typeof(Asap.solve_free), SparseMatrixCSC{Float64,Int}, Matrix{Float64})
    end
    return nothing
end

end # module
