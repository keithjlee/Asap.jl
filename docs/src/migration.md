# Migrating from Asap v0.2.x to v1.0

The complete API mapping for downstream code. The legacy v0.2.x lives on in the registry and the `v0.2.2` tag.

## Types

| v0.2.x | v1.0 |
|---|---|
| `Material(E, G, ρ, ν)` | unchanged (now parametric `Material{T}`; also `Material(E, ρ, ν)` derives G) |
| `Section(A, E, G, Ix, Iy, J, ρ)` | `Section(Material(E, G, ρ, ν), A, Ix, Iy, J)` — material first, density lives on the material |
| `TrussSection(A, E, ρ)` | `Section(Material(E, ρ, ν), A)` — axial-only Section |
| — | **new** `RigiditySection(EA, EIx, EIy, GJ, ρA)` for effective-stiffness workflows (cracked concrete) |
| `Node(pos, fixity)` | unchanged (parametric; position is an `SVector{3}`) |
| `TrussNode(pos, fixity3)` | `Node(pos, vcat(fixity3, trues(3)))` — rotations become *inactive* automatically |
| `Element(n1, n2, sec; release)` | `FrameElement(n1, n2, sec; release, rollangle)` |
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
| `model.S` | `model.cache.K` (free×free only; frozen pattern) — sparse-index workarounds are obsolete: use `ModelState`/`Asap.solve` |
| `model.freeDOFs` etc. | `model.cache.partition.free/.fixed/.inactive` |
| `element.R`, `element.LCS` | `local_frame(el)` (rows = local x/y/z axes); nothing cached on elements |
| `element.length` | `length(el)` |
| `element.Q` | condensed local FEFs in `model.cache.q_local[el.index]` (per segment) |

**DOF-space change for truss models**: legacy truss models used 3 DOFs/node; v1.0 always numbers 6 (inactive rotations never enter the solve). Index a node's translations at `6(i−1)+1 : 6(i−1)+3` — e.g. legacy `u[2:3:end]` (vertical of a 3-DOF layout) becomes `u[2:6:end]`... check each site: the robust form is `displacement(res, node)[2]`.

## Loads

| v0.2.x | v1.0 |
|---|---|
| `NodeForce(node, v)` / `NodeMoment` | unchanged (+ `id`/`case` keywords) |
| `LineLoad(el, v)` | unchanged (lowers to `DistributedLoad`) |
| `PointLoad(el, t, v)` | unchanged. NOTE documented deviation: the axial component now distributes by lever rule, not 50/50 |
| `GravityLoad(el, factor)` | `SelfWeight(el; g, factor)` — driven by `ρA(section)`; the legacy version was broken |
| — | **new** `TrapezoidLoad`, `DistributedLoad` (piecewise-linear, partial-span), `PointMoment`; `case = :tag` on everything |

## Analysis

- `process!` / `solve!(model; reprocess)` unchanged in spirit. Repeated solves reuse the frozen sparsity pattern; `reprocess = true` only after topology changes.
- `solve!(model, loads)` (fresh-load overload) → replace with load `case` tags + `solve_cases!`/`combine` (one factorization, superposition).
- **new** `solve_cases!`, `LoadCombination`, `combine`, `envelope`.
- **new** differentiable path: `state = extract_state(model)` (or build one from design variables), `u = Asap.solve(model, state)`, `compliance(model, state)`. Loading Zygote activates the rule extension; gradients w.r.t. positions, section properties, and connection stiffnesses flow with no further code.

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

Scalar queries with no sampling: `moment_z(model, el, t)`, `shear_y(...)`, etc. (`t` a fraction of length; works uniformly on `VariableElement`). Deflected shapes: `local_displacements(internal_forces(model, el), t)`.

## AsapToolkit absorption

Most of AsapToolkit moved INTO Asap — `using Asap` now provides:

- All parametric generators (`Warren2D`, `Pratt2D`, `BakerTruss`, `TrussFrame`, `SpaceFrame`, `SpaceFrameBeam`, `Frame`, `GridFrame`, `GridNetwork`, ground structures + `to_truss`/`to_frame`). The variable-depth `SpaceFrame` now takes ANY callable `surface(u, v)` on `[0,1]²` — no Interpolations.jl requirement (but its interpolators still work as-is).
- `Geo`/`ModelGeo`/`TrussGeo`/`NetworkGeo`, `ElementDisplacements`, `displacements`.
- `to_network(model)` (FDM translation).
- `clear_supports!`, `element_connectivity`.

AsapToolkit retains only: `SteelSections` (AISC database, `toASAPframe`/`toASAPtruss`), `AsapSections` (polygon geometry), and IO (`topologize`, `GHsave`). There are no re-exports of the moved names — scripts that used them via AsapToolkit should load Asap (most already do).

Other porting notes:

- Release-type dicts (`Asap.FixedFixed` etc. as keys): release types no longer exist as types — key on `release_symbol(el.ends)` instead.
- `element.LCS[i]` → `local_frame(el)[i, :]`; `Asap.nodeids(el)` → `(el.nodeStart.index, el.nodeEnd.index)`.
