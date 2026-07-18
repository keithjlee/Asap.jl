# Solver backends

```@setup solvers
using Asap
steel  = Material(200e6, 77e6, 8.0, 0.3)
wshape = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)
sec = Section(Steel_kNm, 1e-3, 1e-6, 1e-6, 1e-6)
model = Warren2D(11, 1.5, 2.0, sec; load = [0.0, -20.0, 0.0]).model
```

The built-in solver (CHOLMOD Cholesky with LDLᵀ fallback) needs nothing and is the default — most users never think about this. Loading [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) unlocks its entire algorithm collection through one keyword:

```@example solvers
using LinearSolve

solve!(model; solver = KLUFactorization())      # alternative direct solver
solve!(model)                                   # choice is remembered on the model
model.cache.factorization.solver
```

Iterative solvers suit large models:

```@example solvers
p1 = Node([0.0, 0.0, 0.0], :fixed)
p2 = Node([0.0, 0.0, 3.0], :free)
pm = Model([p1, p2], [FrameElement(p1, p2, wshape)],
    [NodeForce(p2, [10.0, 0.0, 0.0])])
solve!(pm; solver = KrylovJL_CG())              # iterative — large models
displacement(pm.results, p2)
```

Repeated solves reuse the factorization's symbolic analysis (numeric-only refactorization on the frozen sparsity pattern) on every backend, including the default. Two notes: unpreconditioned iterative solvers want reasonably conditioned systems — supply a preconditioner for large/stiff models; and on the differentiable path, iterative solvers make gradients inexact adjoints (accuracy follows the solve tolerance) while direct factorizations stay exact.

(`solve!`/`solve` extend the CommonSolve verbs, so loading LinearSolve or other SciML packages never shadows Asap's API.)
