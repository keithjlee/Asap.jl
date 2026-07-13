# Migrating from Asap v0.2.x to v1.0

The complete API mapping for downstream code. The v1.0 core is on branch
`modernization/phase-1`; AsapOptim's port (branch `asap-v1`) is a worked
example of a full migration.

## Types

| v0.2.x | v1.0 |
|---|---|
| `Material(E, G, ρ, ν)` | unchanged (now parametric `Material{T}`; also `Material(E, ρ, ν)` derives G) |
| `Section(A, E, G, Ix, Iy, J, ρ)` | `Section(Material(E, G, ρ, ν), A, Ix, Iy, J)` — material first, density lives on the material |
| `TrussSection(A, E, ρ)` | `Section(Material(E, ρ, ν), A)` — axial-only Section |
| — | **new** `RigiditySection(EA, EIx, EIy, GJ, ρA)` for effective-stiffness workflows (cracked concrete) |
| `Node(pos, fixity)` | unchanged (parametric; position is an `SVector{3}`) |
| `TrussNode(pos, fixity3)` | `Node(pos, vcat(fixity3, trues(3)))` — rotations become *inactive* automatically |
| `Element(n1, n2, sec; release)` | `FrameElement(n1, n2, sec; release, Ψ)` |
| `TrussElement(n1, n2, sec)` | unchanged |
| `BridgeElement` | **deleted** — use `VariableElement(n1, n2, sections, breaks)` (interior joints are internal DOFs, the model is never mutated) |
| `Model(nodes, els, loads)` / `TrussModel(...)` | one `Model` for everything — truss & frame mix freely |
| — | **new** `NodalSpring(node, k6)` elastic supports; `EndSprings`/`EndConditions` semi-rigid joints |

## Results moved OFF the structs

This is the one systematic change every consumer hits:

| v0.2.x | v1.0 |
|---|---|
| `model.u` | `model.results.u` (still full-space; 6 slots/node ALWAYS — truss models included) |
| `model.reactions` | `model.results.reactions` |
| `model.compliance` | `model.results.compliance` |
| `node.displacement` | `displacement(model.results, node) -> SVector{6}` |
| `node.reaction` | `reaction(model.results, node) -> SVector{6}` |
| `element.forces` | `element_forces(model.results, el)` — ALWAYS a 12-vector (truss too) |
| truss axial `forces[2]` | `axial_force(model.results, el)` (slot 7 uniformly) |
| `model.S` | `model.cache.K` (free×free only; frozen pattern) — AsapOptim-style `all_inz` hacks are obsolete: use `ModelState`/`Asap.solve` |
| `model.freeDOFs` etc. | `model.cache.partition.free/.fixed/.inactive` |
| `element.R`, `element.LCS` | `local_frame(el)` (rows = local x/y/z axes); nothing cached on elements |
| `element.length` | `length(el)` |
| `element.Q` | condensed local FEFs in `model.cache.q_local[el.index]` (per segment) |

**DOF-space change for truss models**: legacy truss models used 3 DOFs/node;
v1.0 always numbers 6 (inactive rotations never enter the solve). Index a
node's translations at `6(i−1)+1 : 6(i−1)+3` — e.g. legacy `u[2:3:end]`
(vertical of a 3-DOF layout) becomes `u[2:6:end]`... check each site: the
robust form is `displacement(res, node)[2]`.

## Loads

| v0.2.x | v1.0 |
|---|---|
| `NodeForce(node, v)` / `NodeMoment` | unchanged (+ `id`/`case` keywords) |
| `LineLoad(el, v)` | unchanged (lowers to `DistributedLoad`) |
| `PointLoad(el, t, v)` | unchanged. NOTE documented deviation: the axial component now distributes by lever rule, not 50/50 |
| `GravityLoad(el, factor)` | `SelfWeight(el; g, factor)` — driven by `ρA(section)`; the legacy version was broken |
| — | **new** `TrapezoidLoad`, `DistributedLoad` (piecewise-linear, partial-span), `PointMoment`; `case = :tag` on everything |

## Analysis

- `process!` / `solve!(model; reprocess)` unchanged in spirit. Repeated solves
  reuse the frozen sparsity pattern; `reprocess = true` only after topology
  changes.
- `solve!(model, loads)` (fresh-load overload) → replace with load `case`
  tags + `solve_cases!`/`combine` (one factorization, superposition).
- **new** `solve_cases!`, `LoadCombination`, `combine`, `envelope`.
- **new** differentiable path: `state = extract_state(model)` (or build one
  from design variables), `u = Asap.solve(model, state)`, `compliance(model, state)`.
  Loading Zygote activates the rule extension; gradients w.r.t. positions,
  section properties, and connection stiffnesses flow with no further code.

## Internal forces (was AsapToolkit's ForceAnalysis)

`InternalForces(element, model)` moves INTO core with **axis-correct names**:

| Toolkit v0.1.x | Asap v1.0 |
|---|---|
| `InternalForces(el, model; resolution).x` | `InternalForces(model, el; resolution).x` (note argument order) |
| `.P` | `.N` (tension +) |
| `.My` (was actually about local z!) | `.Mz` |
| `.Vy` | `.Vy` |
| `.Mz` (was about local y) | `−1 × .My` (sign convention verified against pinned fixtures) |
| `.Vz` | `.Vz` |
| — | `.Mx` (torsion — new) |
| `load_envelopes` | `envelope(model, el, caseresults, combos)` — superposition, no re-solves |

Scalar queries with no sampling: `moment_z(model, el, t)`, `shear_y(...)`,
etc. (`t` a fraction of length; works uniformly on `VariableElement`).
Deflected shapes: `local_displacements(internal_forces(model, el), t)`.

## Per-package notes

### AsapToolkit (Phase 5b)
- DELETE `src/ForceAnalysis/` — consumers switch per the table above.
- Direct field reads to update: `.forces` → `element_forces(res, el)`;
  `node.displacement`/`.reaction` → accessors; `element.LCS[i]` →
  `local_frame(el)[i, :]`; `Asap.nodeids(el)` → `(el.nodeStart.index, el.nodeEnd.index)`.
- `Generation/` builders: swap `TrussNode`/`TrussSection`/`Model` constructors
  per the type table; generators returning `model.nodes[:id]` symbol indexing
  still works.
- Release-type dicts (`Asap.FixedFixed` etc. as keys): release types no longer
  exist as types — key on `release_symbol(el.ends)` instead.
- VERIFY: regenerate a few plots side-by-side (visual check).

### AsapHarmonics (Phase 5c)
- `node.reaction` → `reaction(model.results, node)`; `TrussNode`/`TrussModel`
  per table; `e.LCS[1]` → `local_frame(e)[1, :]`; `axial_force(model, el)` →
  `axial_force(model.results, el)`.

### DemandTransport2 (Phase 5d)
- Struct fields typed `::TrussModel` → `::Model{Float64}`; `ElementDemand`'s
  `elements[i].forces[2]` → `axial_force(model.results, elements[i])`.
- AsapOptim consumers: `TrussOptParams` still works (alias of `OptParams`);
  `solve_truss(x, p)` unchanged; `res.U` is now 6 DOFs/node (see above);
  `res.L`, `axial_force(res, p)`, `axial_stress(res, p)`,
  `GeometricProperties(x, p).L/.A`, `updatemodel(p, x)` all preserved.
- The FDM network optimization path (`solve_network`, `QVariable`,
  `NetworkOptParams`) is NOT yet ported — `spaceframe_to_vault` stays on the
  v0.1.x AsapOptim until it lands.
