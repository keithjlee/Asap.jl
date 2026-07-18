# Asap.jl

![](assets/forces-axo.png)

Asap is...
- Another Structural Analysis Package
- results As Soon As Possible
- Analysis of Structures Avec Programming
- A Simple Analysis Please

**Fast, differentiable structural analysis for trusses and frames in Julia** — designed first-and-foremost for information-rich data structures and ease of querying, but always with performance in mind.

Asap solves linear elastic truss and frame structures with the direct stiffness method, recovers exact internal force and displacement fields along members, handles load cases and combinations from a single factorization, and exposes a pure functional analysis path that automatic differentiation engines traverse natively — the same engine that powers gradient-based structural optimization through [AsapOptim](https://github.com/keithjlee/AsapOptim).

See also: [AsapToolkit](https://github.com/keithjlee/AsapToolkit), [AsapHarmonics](https://github.com/keithjlee/AsapHarmonics).

## Installation

```julia
pkg> add Asap
```

## Quick start

```@example quickstart
using Asap

steel = Material(200e6, 80.0, 0.3)                # E [kN/m²], ρ [kg/m³], ν
column = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)   # A, Ix, Iy, J

n1 = Node([0.0, 0.0, 0.0], :fixed)
n2 = Node([0.0, 0.0, 3.0], :free)
el = FrameElement(n1, n2, column)

model = Model([n1, n2], [el], [NodeForce(n2, [10.0, 0.0, 0.0])])
solve!(model)

displacement(model.results, n2)    # SVector{6}: (ux, uy, uz, θx, θy, θz)
```

```@example quickstart
moment_z(model, el, 0.0)           # bending moment at the base
```

## Design

Three separated layers, so each does one job:

1. **Definition** — `Model`, `Node`, elements, loads, springs: what you build. Pure data; nothing analysis-related is cached on it.
2. **Analysis structure** — built once per topology by `process!`: DOF classification, a *frozen* sparsity pattern, scatter maps, a cached factorization. Repeated solves reuse everything.
3. **Results** — returned by `solve!`, queried through accessor functions. Results never live on nodes or elements — which is also what lets dual numbers flow through the differentiable path.

Asap is **units-agnostic**: pick any consistent system (kN–m, N–mm, kip–in) and use it everywhere. Docstrings state dimensions in bracket notation (`[force/length²]`).

Every type prints an engineer-readable summary in the REPL — try typing `steel` or `el` from the example above.

## Where to go from here

The **Guide** walks through every feature with executable examples — start with [Materials and sections](guide/materials-sections.md) and read in order, or jump straight to the topic you need. Coming from the legacy v0.2.x API? See [Migrating from v0.2](migration.md). The **API reference** documents every exported type and function, grouped by concept.

## Verification

The test suite (2,500+ tests) pins: textbook solutions (Kassimali), the exact numerics of the legacy v0.2.x implementation (characterization fixtures), the DiffAnalysis_2024 publication structures and gradients, classic closed-form beam formulas, calculus identities of the recovered force fields, AD-vs-finite-difference gradient checks, and zero-allocation guarantees on the hot paths.

## Citing

When using or extending this software for research purposes, please cite using the following:

```
@software{lee_2024_10581560,
  author       = {Lee, Keith Janghyun},
  title        = {Asap.jl},
  month        = jan,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v0.1},
  doi          = {10.5281/zenodo.10581560},
  url          = {https://doi.org/10.5281/zenodo.10581560}
}
```

Or find a pre-written citation in the style of your choice [here](https://zenodo.org/records/10724610) (see the Citation box on the right side). E.g., for APA:

```
Lee, K. J. (2024). Asap.jl (v0.1). Zenodo. https://doi.org/10.5281/zenodo.10581560
```

## Module reference

```@docs
Asap
```
