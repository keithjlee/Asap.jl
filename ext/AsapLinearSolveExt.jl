"""
    AsapLinearSolveExt

Solver-seam backend for [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl):
loading LinearSolve lets any of its algorithms drive Asap's linear solves,

    using Asap, LinearSolve
    solve!(model; solver = KLUFactorization())
    solve!(model; solver = KrylovJL_CG(precs = ...))   # iterative + preconditioner

The algorithm's `init` cache is held in Asap's `FactorizationCache`, so
repeated solves reuse symbolic analysis / preconditioners the same way the
built-in CHOLMOD path reuses its factorization. Notes:

- The stiffness matrix is passed as-is (numerically symmetric positive
  (semi-)definite on the free DOFs); pick algorithms accordingly.
- Iterative algorithms make gradients through the differentiable path
  *inexact adjoints* — accuracy follows the solve tolerance. Direct
  factorizations remain exact.
- Generic element types (`Model{BigFloat}`, duals) can use LinearSolve's
  generic direct solvers, which the built-in CHOLMOD path cannot.
"""
module AsapLinearSolveExt

using Asap, LinearSolve, SparseArrays, LinearAlgebra

const _Alg = LinearSolve.SciMLLinearSolveAlgorithm

function Asap._factorize(alg::_Alg, K::SparseMatrixCSC)
    prob = LinearProblem(K, zeros(eltype(K), size(K, 1)))
    lincache = init(prob, alg)
    return Asap.FactorizationCache(alg, lincache, nothing)
end

function Asap._refactorize!(alg::_Alg, fc::Asap.FactorizationCache, K::SparseMatrixCSC)
    fc.F.A = K            # marks the LinearSolve cache stale → the next
    return fc             # solve re-runs the numeric factorization
end

function Asap._backsolve(alg::_Alg, fc::Asap.FactorizationCache, b::AbstractVector)
    fc.F.b = b
    sol = LinearSolve.solve!(fc.F)
    return copy(sol.u)    # sol.u aliases the cache's buffer
end

Asap._backsolve(alg::_Alg, fc::Asap.FactorizationCache, B::AbstractMatrix) =
    reduce(hcat, [Asap._backsolve(alg, fc, Vector(col)) for col in eachcol(B)])

end # module
