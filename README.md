[![DOI](https://zenodo.org/badge/426740094.svg)](https://zenodo.org/doi/10.5281/zenodo.10581559)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://keithjlee.github.io/Asap/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://keithjlee.github.io/Asap/dev/)

![](READMEassets/forces-axo.png)

# Asap.jl

Asap is...
- Another Structural Analysis Package
- results As Soon As Possible
- Analysis of Structures Avec Programming
- A Simple Analysis Please

**Fast, differentiable structural analysis for trusses and frames in Julia** — designed first-and-foremost for information-rich data structures and ease of querying, but always with performance in mind.

Asap solves linear elastic truss and frame structures with the direct stiffness method, recovers exact internal force and displacement fields along members, handles load cases and combinations from a single factorization, and exposes a pure functional analysis path that automatic differentiation engines traverse natively — the same engine that powers gradient-based structural optimization through [AsapOptim](https://github.com/keithjlee/AsapOptim).

See also: [AsapToolkit](https://github.com/keithjlee/AsapToolkit), [AsapHarmonics](https://github.com/keithjlee/AsapHarmonics).

# Installation

```julia
pkg> add Asap
```

# Quick start

```julia
using Asap

steel = Material(200e6, 80.0, 0.3)                # E [kN/m²], ρ [kg/m³], ν
column = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)   # A, Ix, Iy, J

n1 = Node([0.0, 0.0, 0.0], :fixed)
n2 = Node([0.0, 0.0, 3.0], :free)
el = FrameElement(n1, n2, column)

model = Model([n1, n2], [el], [NodeForce(n2, [10.0, 0.0, 0.0])])
solve!(model)

displacement(model.results, n2)    # SVector{6}: (ux, uy, uz, θx, θy, θz)
moment_z(model, el, 0.0)           # bending moment at the base
```

# Design

Three separated layers, so each does one job:

1. **Definition** — `Model`, `Node`, elements, loads, springs: what you build. Pure data; nothing analysis-related is cached on it.
2. **Analysis structure** — built once per topology by `process!`: DOF classification, a *frozen* sparsity pattern, scatter maps, a cached factorization. Repeated solves reuse everything.
3. **Results** — returned by `solve!`, queried through accessor functions. Results never live on nodes or elements — which is also what lets dual numbers flow through the differentiable path.

Asap is **units-agnostic**: pick any consistent system (kN–m, N–mm, kip–in) and use it everywhere. Docstrings state dimensions in bracket notation (`[force/length²]`).

Every type prints an engineer-readable summary in the REPL — try typing `steel` or `el` from the example above.

# Features, by example

## Materials and sections

A `Material` carries `E` (Young's modulus), `G` (shear modulus), `ρ` (density), and `ν` (Poisson's ratio):

```julia
steel  = Material(200e6, 77e6, 8.0, 0.3)   # explicit G
timber = Material(13.1e6, 560.0, 0.35)     # isotropic: G derived as E / 2(1+ν)
```

Sections implement a five-accessor contract — `EA`, `EIx`, `EIy`, `GJ`, `ρA` — and the analysis never reads anything else. Two implementations:

**`Section`** — geometry plus material, the standard workflow:

```julia
wshape = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)   # A, Ix (strong), Iy (weak), J
bar    = Section(steel, 5e-3)                     # axial-only: for truss members
EA(wshape), EIx(wshape)                           # derived rigidity products
```

Mismatches are caught early: assigning an axial-only section to a frame element that has moment connections raises a descriptive error at `process!` time, instead of a singular factorization (or a member that silently carries no bending).

**`RigiditySection`** — stores the rigidity *products* directly. This is the natural type for cracked/effective concrete stiffness: linear analysis only ever consumes `EA`, `EIx`, `EIy`, `GJ`, so an equivalent-stiffness section loses nothing.

```julia
Ec, Ig, Ag = 30e6, 8e-4, 0.12
cracked_beam = RigiditySection(
    Ec * Ag,                          # EA
    0.35 * Ec * Ig, 0.35 * Ec * Ig,   # EIx, EIy — ACI 318 cracked-beam modifier
    0.10 * Ec * Ig,                   # GJ
    2.4 * Ag)                         # ρA: mass per length, stored explicitly
```

Any element accepts either type interchangeably.

## Nodes and supports

Every node has six DOFs — three translations, three rotations — with a fixity per DOF (`true` = free). Common supports have symbols; anything else is an explicit 6-vector:

```julia
base   = Node([0.0, 0.0, 0.0], :fixed)      # clamped
pin    = Node([5.0, 0.0, 0.0], :pinned)     # translations fixed, rotations free
roller = Node([10.0, 0.0, 0.0], [true, true, false, true, true, true])  # z-roller
fixnode!(roller, :zfixed)                   # change a support later
```

`planarize!(model)` constrains a model to 2D behavior (fixes out-of-plane DOFs and zeroes roll angles for the `:XY` plane).

Rotational DOFs that nothing stiffens — a node connected only to truss members, or only to fully released member ends — are classified **inactive** and excluded from the solve automatically. No singular stiffness matrices, no manual bookkeeping: an all-truss model solves a translations-only system by construction.

## Frame elements

3D Euler-Bernoulli beam-columns carrying axial force, biaxial bending, shear, and torsion. `rollangle` is the section roll angle about the member axis (default `π/2`).

```julia
girder = FrameElement(n1, n2, wshape, :girder)
rolled = FrameElement(n1, n2, wshape; rollangle = 0.0)
```

### End releases and semi-rigid connections

Connections are **end springs** in the element's local axes. The classical releases are exact limits (`Inf` = rigid, `0` = released), and any finite value is a semi-rigid connection:

```julia
hinged = FrameElement(n1, n2, wshape; release = :fixedfree)   # hinge at far end
# available: :fixedfixed (default), :fixedfree, :freefixed, :freefree, :joist

# semi-rigid: e.g. a bolted end plate with finite rotational stiffness
semi = FrameElement(n1, n2, wshape,
    EndConditions(EndSprings(Inf, Inf, 5e4, 5e4), rigid_end()), :connection)
```

Connection stiffness is ordinary data — it can be a design variable in gradient-based optimization, and everything downstream (fixed-end forces, force recovery) handles it with no special cases.

## Truss elements — and mixing them with frames

Axial-only two-force members. They coexist freely with frame elements in one model; a truss element simply never touches rotational DOFs:

```julia
beam = FrameElement(n1, n2, wshape, :beam)
tie  = TrussElement(n2, n3, bar, :tie)          # n3: an anchor above

model = Model([n1, n2, n3],
    AbstractElement{Float64}[beam, tie],        # mixed vector: type it explicitly
    AbstractLoad{Float64}[NodeForce(n2, [0.0, 0.0, -100.0])])
solve!(model)
axial_force(model.results, tie)                 # tension-positive
```

## Variable elements (varying cross-section)

One user-facing member whose section varies along its length as a chain of prismatic segments — for haunched girders, stepped columns, or optimization results. Interior joints become internal DOFs: the model is never mutated, no phantom nodes appear, and you query the member as a single piece:

```julia
haunched = VariableElement(n1, n2,
    AbstractSection{Float64}[deep, mid, shallow],   # one section per segment
    [0.25, 0.6])                                    # interior break fractions

solve!(model)
moment_z(model, haunched, 0.5)                  # resolves to the right segment
element_forces(model.results, haunched, 2)      # or reach a specific segment
```

`SelfWeight` on a `VariableElement` automatically varies with each segment's `ρA`.

## Spring supports

Elastic supports are **applicative**: a `NodalSpring` references its node (the way a load does) and lives on the model — nodes don't know about their springs, and several springs on one node add up.

```julia
soil = NodalSpring(base, [0.0, 0.0, 5e4, 0.0, 0.0, 0.0], :soil)   # vertical only
pad  = NodalSpring(base, 1e5)                                     # uniform translational
model = Model(nodes, elements, loads; springs = [soil])
```

Spring reactions are recovered as `−k·u` in post-processing.

## Loads

All loads take an `id` and a `case` tag (default `:LC1`):

```julia
NodeForce(n2, [0.0, 0.0, -50.0]; case = :live)
NodeMoment(n2, [0.0, 1e3, 0.0])

LineLoad(beam, [0.0, 0.0, -2.0])                            # uniform, full span
TrapezoidLoad(beam, 0.2, 0.8, 1.0, 3.0, [0.0, 0.0, -1.0])   # partial-span, varying
DistributedLoad(beam, [0.0, 0.3, 0.5], [0.0, 4.0, 0.0],     # arbitrary piecewise-linear
    [0.0, 0.0, -1.0])
PointLoad(beam, 0.4, [0.0, 0.0, -10.0])                     # at 40% of the span
PointMoment(beam, 0.5, [0.0, 0.0, 800.0])                   # concentrated moment

SelfWeight(beam; g = [0.0, 0.0, -9.81])                     # from ρA(section)
```

Every distributed shape lowers to one canonical piecewise-linear type with one exact integration (3-point Gauss against the element shape functions — exact, not approximate, for these loads). Adding a custom element load type means implementing a single method, `fixed_end_forces`.

## Analysis and results

```julia
model = Model(nodes, elements, loads)
solve!(model)                     # process! runs automatically the first time

res = model.results
displacement(res, node)           # SVector{6} at a node
reaction(res, node)               # support reactions (incl. moment reactions)
element_forces(res, el)           # local end-force 12-vector
axial_force(res, el)              # tension-positive scalar
res.compliance                    # external work Fᵀu
```

Repeated solves — geometry or section *values* changed, topology unchanged — reuse the frozen sparsity pattern and refactorize numerically: just call `solve!(model)` again. After topology changes (elements added, releases toggled, supports changed): `solve!(model; reprocess = true)`.

## Internal forces and deflected shapes

Recovery is **equilibrium-based**: starting from the exact member end actions, the applied loads are integrated analytically along the member — the recovered fields are exact for any end condition (releases, semi-rigid springs) with no per-case formulas.

```julia
# scalar queries at any fraction of the length (zero allocation):
moment_z(model, beam, 0.5)     # bending about local z — pairs with Vy, sagging+ for +y loads
shear_y(model, beam, 0.25)
axial_force(model, beam, 0.0)
torsion(model, beam, 0.5)

# dense sampling for plotting — stations include every load breakpoint and
# BOTH sides of each point action, so shear jumps render as true jumps:
f = InternalForces(model, beam; resolution = 40)
f.x; f.N; f.Vy; f.Mz; f.Vz; f.My; f.Mx

# deflected shape: exact local displacements at any fraction
st = internal_forces(model, beam)
local_displacements(st, 0.5)   # (u, v, w) in local axes
```

Names are axis-correct: `Mz` is the moment about local z and pairs with `Vy` (`dMz/dx = Vy`); `My` pairs with `Vz` (`dMy/dx = −Vz`).

## Load cases, combinations, envelopes

Tag loads with cases; solve every case against **one** factorization; get any factored combination by superposition — no re-solves, ever:

```julia
loads = AbstractLoad{Float64}[
    LineLoad(beam, [0.0, 0.0, -2.0]; case = :dead),
    PointLoad(beam, 0.4, [0.0, 0.0, -10.0]; case = :live),
    NodeForce(n2, [5.0, 0.0, 0.0]; case = :wind),
]
model = Model(nodes, elements, loads)

cr = solve_cases!(model)                    # one assembly, one factorization
strength = LoadCombination(:LRFD, [:dead => 1.2, :live => 1.6, :wind => 0.5])
res = combine(cr, strength)                 # exact, by superposition
displacement(res, n2)

# station-wise min/max over a combination set — what design checks consume:
env = envelope(model, beam, cr,
    [strength, LoadCombination(:service, [:dead => 1.0, :live => 1.0])])
env.x; env.lo; env.hi                       # rows: N, Vy, Mz, Vz, My, Mx
```

## Differentiable analysis

The entire pipeline has a **pure functional path**: extract a `ModelState` (positions and sections as plain data), evaluate, differentiate. Loading Zygote (or anything ChainRules-aware) activates Asap's rule extension automatically; a handful of small rules cover the sparse construction and the linear solve, and everything else differentiates natively.

```julia
using Zygote

solve!(model)                               # builds the analysis cache
state = extract_state(model)

# gradient of compliance w.r.t. every node coordinate — a 3 × n sensitivity field:
g = Zygote.gradient(state.X) do X
    compliance(model, ModelState{Float64}(X, state.sections))
end[1]
```

Gradients flow with respect to node positions, section properties, and semi-rigid connection stiffnesses. For design-variable bookkeeping (areas, coupled symmetric geometry, bounds) and optimization-ready objectives, use [AsapOptim](https://github.com/keithjlee/AsapOptim) — a thin layer over this path, verified against the original published implementation to 13 digits (see its `docs/AD_VERIFICATION_AND_BENCHMARKS.md`).

## Force density method (form-finding)

A self-contained FDM subsystem for cable/tension networks ships in the same package, built on the same three-layer architecture as the frame core: `Network` holds pure definition data, `process!` freezes the force-density pattern once per topology, and repeated solves are a pure scatter of the current `q` values plus a numeric-only refactorization — so form-finding sweeps are cheap by construction. `solve!` updates the node positions in place (the geometry *is* the result); forces and reactions are queried from `network.results`:

```julia
# a 2-cable chain: anchors at the ends, load at the middle
a = FDMnode([0.0, 0.0, 0.0], false)          # false = fully fixed (anchor)
m = FDMnode([1.0, 0.0, 0.0], true, :mid)     # true  = fully free
b = FDMnode([2.0, 0.0, 0.0], false)

cables = [FDMelement(a, m, 2.0), FDMelement(m, b, 2.0)]   # q = 2.0 [force/length]
net = Network([a, m, b], cables, [FDMload(m, [0.0, 0.0, -1.0])])

solve!(net)
m.position                            # ([1.0, 0.0, -0.25]) — the found geometry
member_force(net.results, cables[1])  # cable tension q·L
reaction(net.results, a)              # anchor force, per fixed axis

# form-finding sweep: mutate q and re-solve (numeric refactorization only)
update_q!(net, 4.0)                   # halves the sag
```

Fixity is **per-axis** (`true` = free, like `Node`): prescribe the plan and let the height find equilibrium — natural for vault form-finding. FDM equilibrium is separable per coordinate, so each axis solves against its own free/fixed partition:

```julia
crown = FDMnode([1.0, 1.0, 0.0], Bool[false, false, true])   # plan fixed, z free
```

The `solver` keyword works here too (any LinearSolve algorithm once loaded). `to_network(model)` converts a solved truss model into an equivalent FDM network, and `to_truss(network, section)` converts back. Differentiable force-density optimization (`QVariable`, `solve_network`) lives in AsapOptim.

## Parametric structure generators

A family of one-line generators (absorbed from AsapToolkit) builds, loads, and solves common structural typologies — handy for testing, benchmarking, and optimization studies:

```julia
sec = Section(Steel_kNm, 1e-3, 1e-6, 1e-6, 1e-6)

truss = Warren2D(11, 1.5, 2.0, sec; load = [0.0, -20.0, 0.0])   # 2D Warren arch
grid  = SpaceFrame(6, 1.2, 6, 1.2, 1.0, sec; support = :corner) # double-layer grid

# variable-depth spaceframe: any callable surface(u, v) on [0,1]² sets the top layer
vault = SpaceFrame(6, 1.2, 6, 1.2, 1.0, (u, v) -> 0.5 * sinpi(u) * sinpi(v), sec)

# every generator returns a solved model plus its generation parameters
maximum(abs, truss.model.results.u)

# 3D building frame: columns, primary beams, joists (with releases), braces
fsec = Section(Steel_kNm, 1e-2, 1e-4, 5e-5, 1e-6)
bldg = Frame(2, 6.0, 2, 5.0, 2, 4.0, 2.0, fsec, fsec, fsec, fsec)

# ground structures for layout optimization: dense candidate grids → models
gs = XGroundStructure(6.0, 4, 4.0, 3)
candidates = to_truss(gs, sec; load = [1.0, 0.0, 0.0])
```

The full set: `Warren2D`, `Pratt2D`, `BakerTruss`, `TrussFrame`, `SpaceFrame`, `SpaceFrameBeam`, `Frame`, `GridFrame`, `GridNetwork` (FDM), and the `XGroundStructure` / `DenseGroundStructure` / `BoundedGroundStructure` family with `to_truss` / `to_frame`.

## Plot-ready geometry extraction

`Geo(model)` (or `Geo(network)`) flattens a solved structure into plain arrays for plotting — node positions, displaced positions, element connectivity, and per-element force ranges with their maxima. `ElementDisplacements(element, model)` samples the exact deflected curve along a member:

```julia
geo = Geo(truss.model)          # geo.nodes, geo.disp, geo.indices, geo.P, …

ed = ElementDisplacements(bldg.model.elements[1], bldg.model; resolution = 20)
ed.basepositions .+ 100 .* ed.uglobal   # 3×20 displaced curve, scaled ×100
```

## Solver backends

The built-in solver (CHOLMOD Cholesky with LDLᵀ fallback) needs nothing and is the default — most users never think about this. Loading [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) unlocks its entire algorithm collection through one keyword:

```julia
using LinearSolve

solve!(model; solver = KLUFactorization())      # alternative direct solver
solve!(model; solver = KrylovJL_CG())           # iterative — large models
solve!(model)                                   # choice is remembered on the model
```

Repeated solves reuse the factorization's symbolic analysis (numeric-only refactorization on the frozen sparsity pattern) on every backend, including the default. Two notes: unpreconditioned iterative solvers want reasonably conditioned systems — supply a preconditioner for large/stiff models; and on the differentiable path, iterative solvers make gradients inexact adjoints (accuracy follows the solve tolerance) while direct factorizations stay exact.

(`solve!`/`solve` extend the CommonSolve verbs, so loading LinearSolve or other SciML packages never shadows Asap's API.)

# Documentation map

- **[keithjlee.github.io/Asap](https://keithjlee.github.io/Asap/stable/)** — the full documentation site: guide pages for every feature (with executed examples), the v0.2 → v1.0 migration guide, and the complete API reference.
- `docs/releases/` — release notes, including the v0.2.x → v1.0 breaking-change summary
- Every exported function and type carries a full docstring (`?Section`, `?dof_signature`, …), and every type has an engineer-readable REPL display.

# Verification

The test suite (2,500+ tests) pins: textbook solutions (Kassimali), the exact numerics of the legacy v0.2.x implementation (characterization fixtures), the DiffAnalysis_2024 publication structures and gradients, classic closed-form beam formulas, calculus identities of the recovered force fields, AD-vs-finite-difference gradient checks, and zero-allocation guarantees on the hot paths.

# Citing

When using or extending this software for research purposes, please cite using the following:

### Bibtex
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

### Other styles
Or find a pre-written citation in the style of your choice [here](https://zenodo.org/records/10724610) (see the Citation box on the right side). E.g., for APA:
```
Lee, K. J. (2024). Asap.jl (v0.1). Zenodo. https://doi.org/10.5281/zenodo.10581560
```
