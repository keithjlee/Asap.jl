# Differentiable path

The pure functional analysis path: `ModelState`, `extract_state`, `solve`, and `compliance`.

```@autodocs
Modules = [Asap]
Pages = ["Analysis/functional.jl"]
Private = false
```

## Load cases

Several load cases at one design state, one factorization: the plain-data case table and its
design-dependent fixed-end forces.

```@autodocs
Modules = [Asap]
Pages = ["Analysis/case_loads.jl"]
Private = false
```

## Internals

The linear-solve seam of the differentiable path (the function AD rules attach to):

```@docs
Asap.solve_free
```

The plain-data mirror of the applied element loads that [`assemble_Pf`](@ref) re-evaluates per state:

```@docs
Asap.FEFTable
```
