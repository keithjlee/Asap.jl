# Differentiable path

The pure functional analysis path: `ModelState`, `extract_state`, `solve`, and `compliance`.

```@autodocs
Modules = [Asap]
Pages = ["Analysis/functional.jl"]
Private = false
```

## Internals

The linear-solve seam of the differentiable path (the function AD rules attach to):

```@docs
Asap.solve_free
```
