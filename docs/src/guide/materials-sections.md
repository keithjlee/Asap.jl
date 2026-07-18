# Materials and sections

```@setup materials
using Asap
```

A [`Material`](@ref) carries `E` (Young's modulus), `G` (shear modulus), `ρ` (density), and `ν` (Poisson's ratio):

```@example materials
steel  = Material(200e6, 77e6, 8.0, 0.3)   # explicit G
timber = Material(13.1e6, 560.0, 0.35)     # isotropic: G derived as E / 2(1+ν)
```

Sections implement a five-accessor contract — [`EA`](@ref), [`EIx`](@ref), [`EIy`](@ref), [`GJ`](@ref), [`ρA`](@ref) — and the analysis never reads anything else. Two implementations:

## `Section` — geometry plus material

The standard workflow:

```@example materials
wshape = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)   # A, Ix (strong), Iy (weak), J
bar    = Section(steel, 5e-3)                     # axial-only: for truss members
EA(wshape), EIx(wshape)                           # derived rigidity products
```

Mismatches are caught early: assigning an axial-only section to a frame element that has moment connections raises a descriptive error at `process!` time, instead of a singular factorization (or a member that silently carries no bending).

## `RigiditySection` — rigidity products directly

[`RigiditySection`](@ref) stores the rigidity *products* directly. This is the natural type for cracked/effective concrete stiffness: linear analysis only ever consumes `EA`, `EIx`, `EIy`, `GJ`, so an equivalent-stiffness section loses nothing.

```@example materials
Ec, Ig, Ag = 30e6, 8e-4, 0.12
cracked_beam = RigiditySection(
    Ec * Ag,                          # EA
    0.35 * Ec * Ig, 0.35 * Ec * Ig,   # EIx, EIy — ACI 318 cracked-beam modifier
    0.10 * Ec * Ig,                   # GJ
    2.4 * Ag)                         # ρA: mass per length, stored explicitly
```

Any element accepts either type interchangeably.
