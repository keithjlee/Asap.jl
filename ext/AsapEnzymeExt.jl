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
using Enzyme: EnzymeRules, Annotation, Const, Duplicated, BatchDuplicated
using SparseArrays

# Enzyme's @import_rrule is implemented inside ITS ChainRulesCore extension,
# which is not reliably active while THIS extension precompiles — register
# the imported rules at load time instead.
function __init__()
    @eval begin
        Enzyme.@import_rrule(typeof(Asap.solve_free), SparseMatrixCSC{Float64,Int}, Vector{Float64})
        Enzyme.@import_rrule(typeof(Asap.solve_free), SparseMatrixCSC{Float64,Int}, Matrix{Float64})
        Enzyme.@import_rrule(typeof(Asap.sparse_from_pattern), SparseMatrixCSC{Float64,Int}, Vector{Float64})
    end
    return nothing
end

#=
FORWARD mode: a NATIVE EnzymeRules.forward rule for solve_free — the
implicit-function theorem `K u̇ = Ḟ − K̇ u`, one back-substitution per
tangent on the forward factorization.

Native (not `@import_frule`) for a correctness-critical reason: under
runtime activity, Enzyme passes runtime-CONSTANT arguments as Duplicated
with the shadow ALIASING the primal (`dval === val`). The generic frule
bridge cannot know that convention and reads the alias as a real tangent
(`ΔF = F`), producing silently wrong derivatives. This rule checks the
aliasing explicitly and treats those arguments as constants.
=#

# a tangent is real only if the argument is Duplicated AND its shadow is
# distinct storage (runtime-activity aliasing ⇒ runtime constant)
_shadow_or_nothing(arg::Const, _) = nothing
_shadow_or_nothing(arg::Duplicated, _) = arg.dval === arg.val ? nothing : arg.dval
_shadow_or_nothing(arg::BatchDuplicated, i::Int) =
    arg.dval[i] === arg.val ? nothing : arg.dval[i]

EnzymeRules.forward(config::EnzymeRules.FwdConfig,
    ::Const{typeof(Asap.solve_free)}, RT::Type,
    K::Annotation{<:SparseMatrixCSC{Float64,Int}},
    F::Annotation{<:Union{Vector{Float64},Matrix{Float64}}}) =
    _solve_free_fwd(config, RT, nothing, K, F)

# 3-arg (solver seam) — the solver object is opaque configuration, treated
# as data regardless of its annotation (runtime activity hands MUTABLE
# solvers like CachedSolver over as Duplicated with an aliased shadow, so
# a Const-only signature would silently not match and Enzyme would descend
# into the factorization internals)
EnzymeRules.forward(config::EnzymeRules.FwdConfig,
    ::Const{typeof(Asap.solve_free)}, RT::Type, solver::Annotation,
    K::Annotation{<:SparseMatrixCSC{Float64,Int}},
    F::Annotation{<:Union{Vector{Float64},Matrix{Float64}}}) =
    _solve_free_fwd(config, RT, solver.val, K, F)

function _solve_free_fwd(config, RT::Type, solver, K, F)
    Kp = K.val
    fact = Asap._factorize(solver, Kp)
    # FactorizationCache fields are deliberately untyped — assert the solve
    # results so the rule returns concretely-typed (Batch)Duplicated values
    u = (fact \ F.val)::typeof(F.val)

    needs_shadow = EnzymeRules.needs_shadow(config)
    needs_primal = EnzymeRules.needs_primal(config)
    needs_shadow || return needs_primal ? u : nothing

    function tangent(i)
        K̇ = _shadow_or_nothing(K, i)
        Ḟ = _shadow_or_nothing(F, i)
        rhs = Ḟ === nothing ? zero(u) : copy(Ḟ)
        K̇ === nothing || (rhs .-= K̇ * u)
        return (fact \ rhs)::typeof(F.val)
    end

    # Val-typed width: EnzymeRules.width(config) constant-folds from the
    # config type, so the tangent tuple is a concrete NTuple. The returned
    # wrapper must match RT EXACTLY — a type-unstable call site infers
    # RT = (Batch)Duplicated{Any, W}, and Enzyme rejects the concretely
    # typed wrapper even though it subtypes it — so construct through RT
    # whenever RT is fully parameterized.
    W = Val(EnzymeRules.width(config))
    if W === Val(1)
        du = tangent(1)
        needs_primal || return du
        DT = RT isa UnionAll ? Duplicated : RT
        return DT(u, du)
    else
        dus = ntuple(tangent, W)
        needs_primal || return dus
        BT = RT isa UnionAll ? BatchDuplicated : RT
        return BT(u, dus)
    end
end

end # module
