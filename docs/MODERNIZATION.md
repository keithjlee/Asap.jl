# Asap.jl v1.0 Modernization Plan

> Living roadmap — update the phase table as work lands.
> Status: **Phase 0 in progress** (2026-07-13).

## Context

Asap.jl (v0.2.2) is Keith's structural FEA package, 5–6 years old, Float64-strict, with fully duplicated truss/frame code paths, allocation-heavy assembly (COO rebuilt every solve — the primary bottleneck), fragile BridgeElement "shatter" processing, binary releases that leave singular torsional modes, and internal-force diagrams living in AsapToolkit as superposed analytic formulas. AsapOptim (the feature repo) re-implements the whole pipeline functionally with hand-written Zygote rrules and reverse-engineers `model.S`'s CSC internals. Goal: modernize Asap into the fastest open-source FEA engine for structural research (long-term: displace OpenSees), integrate into Forma for full-building concrete analysis, and collapse AsapOptim's workaround layer.

## Locked decisions (Keith, 2026-07-13)

- **Clean break**: new major version; AsapOptim/AsapToolkit/AsapHarmonics/DemandTransport2 updated in lockstep, no shims.
- **Branch policy**: `main` stays the stable v0.2.x for ongoing projects. ALL modernization work happens on `modernization/*` branches (Phase 0 on `modernization/phase-0`; later phases branch from the previous phase). Merge to `main` only when Keith calls v1.0 ready.
- **AD-first core**: Asap provides both an in-place fast path and a pure functional path sharing one element kernel. AsapOptim shrinks to variables/indexers/objectives.
- **Julia 1.10 LTS** floor. Deps: LinearAlgebra, SparseArrays, **+ StaticArrays** (hard), **+ ChainRulesCore** (weak dep / extension). No Unitful. No LinearSolve for now (CHOLMOD `cholesky!` behind a `factorize!` seam; extension later).
- **First feature after core rework**: shape-function/equilibrium internal-force recovery.
- Style stays "Keith-like": readable, explicit, directory-organized.

## Architecture (synthesized from both designs)

### Three-layer separation
1. **Definition** — `Model{T}`: nodes, elements, sections, loads, springs. Pure user data, no analysis state.
2. **Analysis structure** — `AnalysisCache{T}` built once by `process!`: DOF partition, element groups, frozen sparsity pattern, nzmap/scatter, buffers, cached factorization.
3. **Results** — `LinearResults{T}`: full-space u, reactions, element end forces, compliance. Result scalar type independent of model's T (Duals flow through the pure path). Accessors replace struct mutation: `displacement(res, node)`, `element_forces(res, el)`, `moment_at(res, el, t)`.

### Types
- `Material{T}`: E, G, ρ, ν (G derivable from E, ν).
- `AbstractSection{T}` with **accessor contract** `EA/EIx/EIy/GJ/ρA(s)` — kernels never read fields:
  - `Section{T}`: material + A, Ix, Iy, J (products computed in accessors).
  - `RigiditySection{T}`: stores EA, EIx, EIy, GJ, ρA **directly** — cracked concrete (ACI 318 modifiers, Branson), homogenized members, optimization parameterizations. TrussSection deleted (truss only calls EA/ρA).
- `Node{T}`: SVector position, 6-DOF fixity, always 6 slots. **TrussNode deleted** — rotational DOFs become *inactive* via element participation.
- Elements — analysis state (K, Q, R, LCS, length, globalID, forces) **removed from structs**:
  - `FrameElement{T,S}` (rename of Element): nodes, section, `EndConditions{T}`, Ψ.
  - `TrussElement{T,S}`.
  - `VariableElement{T,S}`: one user-facing element = chain of prismatic segments (`sections::Vector{S}`, `breaks`), via **internal DOFs** (6 per interior joint, allocated after nodal DOFs by process!) — NOT model-node shattering, NOT static condensation (keeps mass/geometric-stiffness exact later; recovery = same code path as primitives). Replaces BridgeElement; `bridgeprocessing.jl` deleted.
- **Releases → end-spring data**: `EndSprings{T}` (kx, kt, ky, kz per end; Inf=rigid, 0=released). Release symbols (`:fixedfixed`…`:joist`) map to exact limits and stay the user API. Element vectors become concretely typed (releases no longer type params). Connection stiffness becomes a **differentiable design variable** (Forma concrete workflows).
- `NodalSpring{T}`: **applicative data, not a node property** — a standalone struct referencing a node (like a load references its target), stored in `model.springs`. Diagonal stiffness per global DOF; multiple springs on one node compose additively; marks DOFs active; reaction −k·u. Nodes never know about their springs.
- Loads: immutable, parametric, `case::Symbol` tag, **no release type param**. `NodeForce/NodeMoment`; `DistributedLoad{T}` = canonical piecewise-linear intensity (breakpoints t ∈ [0,1], w per breakpoint, direction, :global/:local) with `LineLoad`/`TrapezoidLoad`/`TributaryLoad` convenience constructors lowering to it; `PointLoad`; new `PointMoment`; `SelfWeight` (lowers to DistributedLoad via `ρA(section)`, per-segment for VariableElement — GravityLoad's ρ bug impossible by construction).

### Element interface (extensibility contract)
`nodes(el)`, `n_internal_dofs(el)`, `ndofs(el)`, `dof_signature(el)` (which DOF slots the element actually couples — drives activity), `local_frame(el, x1, x2)`, `stiffness(el, x1, x2) -> SMatrix` (pure kernel, positions passed explicitly so fast + AD paths share it), `fixed_end_forces`, `local_forces`; future hooks `mass`, `geometric_stiffness` (declared now, same assembly machinery later).

### DOF management (fixes the hinge/torsion singularity structurally)
Global space = 6N nodal slots + internal blocks. Three states per DOF: **free / fixed / inactive**. Activity accumulated from element `dof_signature`s + springs: truss-only model → 3N system automatically; a rotation touched only by released ends → inactive, never enters K — no singular mode, no regularization. Mixed truss/frame just works. `process!` validates loads on inactive DOFs with a named error (today: bare SingularException).

### Assembly
- **Symbolic pass (once)**: group elements by concrete type (type-stable function barriers over heterogeneous `model.elements`), build free×free pattern, freeze colptr/rowval, compute per-element `nzmap` (first-class version of AsapOptim's `all_inz` hack) + Boolean `scatter` matrix for the pure path.
- **In-place path**: `fill!(nzval, 0)`, per-group SMatrix kernel → scatter via nzmap; springs added to diagonal. Zero allocations per solve. Λ transform done blockwise (16 3×3 blocks).
- **Pure path**: `ModelState{T}` (positions + sections extracted from model), same kernels mapped non-mutating, `nzval = scatter * V` (native sparse-matvec adjoint), `SparseMatrixCSC(n, n, colptr, rowval, nzval)`. Parity test: pure ≡ in-place.
- **`ext/AsapChainRulesExt.jl`** — the only 3 rrules needed: sparse ctor (project to nzval), cached-factorization linear solve adjoint, optional nonzeros-gather. This replaces AsapOptim's entire adjoint layer.
- K assembled **free×free only**; reactions recovered element-wise (no full-space K ever formed).

### Semi-rigid joints (Monforton–Wu fixity factors)
`pᵢ = 1/(1 + 3EI/(kᵢL))`; 4×4 bending block per plane with `EI/(L³(4−p₁p₂))` prefactor (reproduces `k_fixedfixed` at p=1,1 and `k_fixedfree` at p=1,0 exactly); axial/torsion end springs as series combinations. FEF correction = generic 2×2 static condensation `q_semi = q_n − K_nb·K_bb⁻¹·q_b` — replaces all hand-written `q(::ElementLoad{R})` methods. Old closed-form matrices in `Elements/K.jl` + `fixedEndForces.jl:114-201` become **test oracles**, then die. AD guard: branch `isinf(k)` explicitly (Duals through Inf → NaN).

### Loads / FEF engine
Single extension point: `fixed_end_forces(load, el) -> SVector{12,T}` clamped-clamped local vector; condensation applied generically afterward. Engine: **3-point Gauss quadrature per linear segment against Hermite/linear shape functions** — integrand quartic → exact, ~40 lines replacing the whole closed-form catalog, AD-transparent. Tasha's kernels ported as test oracles only.

### Solve pipeline
`process!` (once per topology) → `solve!(model)` (assemble in place → `cholesky!` numeric refactorization on frozen pattern, `ldlt!` fallback → u → postprocess → `model.results`); `solve!(model, loads)` skips assembly/factorization (one back-substitution); pure `solve(cache, state, loads)` for AD. Stale-cache detection via topology hash. **Load cases/combos**: multi-RHS solve (U matrix ndof×ncases, one factorization), `LoadCombination` results by superposition, envelopes = extrema over combos per station.

### Internal-force recovery (priority feature; lives in core, not Toolkit)
**Equilibrium recovery**: start from exact local end forces `f = k_local·u_local + q̃`, integrate applied loads analytically along the member (`Mz(x) = f₂x − f₆ + ∫w_y(s)(x−s)ds + Σ point terms`). Exact everywhere for consistent-FEF Euler-Bernoulli, any end condition — no per-release formula zoo. `ElementForceState{T}` (end actions + `LoadTrace` with precomputed cumulative load integrals) → scalar zero-allocation evaluators (`moment_z(s,t)`, `shear_y(s,t)`, …) + `sample(s; n) -> InternalForces{T}` (fields x, N, Vy, Mz, Vz, My, Mx) for plotting; stations include breakpoints and both sides of point loads. VariableElement: composite state, binary search super-fraction → segment. Displacement recovery: Hermite homogeneous part + double-integrated exact curvature (piecewise quintic) — replaces Toolkit's ~300-line D/Theta catalog.
**Toolkit port trap**: current `.My` is actually moment about local z. Map: `old.My → new.Mz`, `old.Mz → new.My`, `P → N` + sign audit ([-1,1,1,1,-1] flips at ForceAnalysis.jl:89). Toolkit's Flocal omits `element.Q` → GravityLoad diagrams silently wrong today.

### Concrete stiffness (Keith's question — answered)
**Yes, equivalent E·I suffices**: linear member FEA consumes material/section data only through EA, EIx, EIy, GJ products (every `k_fixedfixed` entry is a product × geometry factor). Transformed-section vs equivalent-E are numerically identical. `RigiditySection` is the abstraction. Limits: cracking varying along member → VariableElement with segment-wise RigiditySections; axial-flexural interaction → iterate section properties; creep → per-case section sets; Timoshenko later → add GAs products.

## Reconciliations between the two designs
- Loads reference **element/node objects** (core design), not integer indices — matches current style; `case::Symbol` added from loads design. Indices (`index::Int`) assigned by `process!`.
- FEF flow: load kernel returns clamped local 12-vector (loads design) → element layer condenses per `EndConditions` and rotates to global (core design). Composable, both stand.
- Section fields: use core naming `EA, EIx, EIy, GJ, ρA` (Ix/Iy matches existing Asap convention).
- Results: single `LinearResults{T}` type generalized to hold `U::Matrix{T}` (ndof × ncases) + retained factorization (loads design's AnalysisResult features fold in at Phase 4).

## New src/ layout

```
src/
  Asap.jl
  Materials/ {materials.jl, sections.jl}
  Nodes/nodes.jl
  Elements/ {interface.jl, end_conditions.jl, frame.jl, truss.jl, variable.jl,
             kernels/{transformation.jl, stiffness.jl, fixed_end_forces.jl, mass.jl}}
  Springs/nodal_springs.jl
  Loads/loads.jl
  Analysis/ {dofs.jl, symbolic.jl, assemble.jl, functional.jl, factorize.jl, solve.jl}
  Results/ {results.jl, postprocess.jl, queries.jl}
  Model/ {model.jl, utilities.jl}
  FDM/            # unchanged for now
ext/AsapChainRulesExt.jl
```

## Phased roadmap

| Phase | Scope | Exit criterion |
|---|---|---|
| **0** | Hygiene & characterization: CLAUDE.md (Asap + 3 satellites), **commit this plan into the repo as `docs/MODERNIZATION.md`** (living roadmap, updated as phases land), CI (1.10 + release), **characterization tests pinning current numerics** (Kassimali examples; FEF `q()` per load×release; element.forces; reactions; u), Toolkit InternalForces fixtures (canonical 3D portal frame, all releases), failing GravityLoad-bug test | Old Asap passes pinned suite; fixtures committed |
| **1** | Core rework: Materials/Sections + accessors → Node{T} → EndConditions (tested vs old K.jl) → element kernels (pure SMatrix) → DOF layer (activity/partition/internal blocks) → symbolic + in-place assembly → solve pipeline → LinearResults → pure path + ext rrules (parity + Zygote-vs-FiniteDiff tests) → VariableElement (delete bridgeprocessing) | Characterization suite passes on new core; parity fast≡pure; gradcheck green |
| **2** | Loads: Gauss/Hermite FEF engine, new taxonomy, generic condensation | New FEFs match pinned q() to 1e-12; trapezoid vs Tasha's kernels; SelfWeight passes GravityLoad-bug test |
| **3** | **Force recovery** (priority): LoadTrace, ElementForceState, evaluators, sample, displacement recovery, VariableElement dispatch | Matches Toolkit fixtures under My/Mz rename; property tests Mz′=Vy, Vy′=−w; `@allocated == 0` scalar evals |
| **4** | Cases/combos/envelopes: multi-RHS solve, LoadCombination, envelope | Envelope ≡ brute-force per-combo solve; single factorization |
| **5a** | AsapOptim: delete Functions/{K,Kframe,Ktruss,Rframe,Rtruss}.jl + assembly half of Solve.jl + all_inz layer; keep Variables/Indexers/Constraints; call `extract_state`/`solve(cache, state)` | Suite green; gradient check vs FiniteDifferences (truss + frame) |
| **5b** | AsapToolkit: delete ForceAnalysis/ (absorbed); mechanical rename port (P→N, My↔Mz, accessor migration); keep Generation/SteelSections/AsapSections/IO | Plots regenerate identically vs fixtures |
| **5c** | AsapHarmonics: rename pass (`node.reaction` → `reaction(res, node)` etc.) | Suite green |
| **5d** | DemandTransport2: port legacy Truss types out of struct fields (`model_opt::TrussModel`, `ElementDemand(::TrussModel)`, `ResponsiveGeo2D(::TrussModel)` → unified `Model`); replace the one struct-embedded result read `element.forces[2]` (`src/Demand/data_structures.jl:47`) with the results accessor; objective functions move to the new differentiable results API | Truss + network optimization examples reproduce prior optima; Zygote gradients green |
| Later | Geometric nonlinearity (`geometric_stiffness` hook + tangent assembly reuse), dynamics (`mass` hook), FDM parametrization, Forma integration | — |

### Downstream AD contract (driven by DemandTransport2, the reference consumer)
DemandTransport2 (`../DemandTransport2`, optimal transport + differentiable structural analysis) is a pure consumer of AsapOptim's differentiable API — no rrules or sparsity hacks of its own. The Phase 5a AsapOptim rework must preserve, on top of the new core pure path:
- Gradients of displacements, axial forces, and element lengths **w.r.t. node positions** (`SpatialVariable`) and of the FDM network solve w.r.t. **force densities** (`QVariable`). It does not optimize areas.
- A results object exposing displacements + element lengths (today `res.U`, `res.L`) and `axial_force(res, p)`; network results `res.X/Y/Z/Q`.
- A **geometry-only, no-solve path** (today `GeometricProperties(x, p)` → `.L`) for length-based objectives.
- `updatemodel(params, x)` returning a materialized, solved model for post-processing/viz.
- Fixed sparsity across repeated solves of a fixed-topology model (already guaranteed by `AnalysisCache`).
Note DemandTransport2 also exercises the **FDM/Network path** (`solve_network`, `NetworkOptParams`) — FDM stays functional through the transition even though its parametrization is deferred. Its own migration is Phase 5d; hardest hits are the legacy Truss types hard-coded in its struct fields and dispatch signatures (fails at precompile, not just call sites).

## Critical files
- `src/Elements/K.jl` — closed-form matrices → generalized fixity-factor kernel (and test oracle)
- `src/Model/preprocessing.jl` — DOF numbering + COO assembly → activity/symbolic/nzmap design
- `src/Model/analysis.jl` — process!/solve! → AnalysisCache + cached factorization
- `src/Loads/fixedEndForces.jl` — release-specific FEFs → Gauss engine + generic condensation (broken GravityLoad dies here)
- `src/Elements/elements.jl` — element family + VariableElement (BridgeElement deleted)
- `src/Model/postprocessing.jl` — seed of ElementForceState.f0
- AsapOptim `Types/Utilities.jl` (all_inz), `Functions/*` — deletion targets in 5a
- AsapToolkit `src/ForceAnalysis/*` — API to preserve (axis-fixed), then absorb

## Known bugs to fix/regress along the way
GravityLoad FEF (nonexistent `load.element.ρ`); dangling `release!` export; `copy(::TrussModel)` (`element.nodeIDs`); `shatterReleaseDict` duplicate key; stale `process_bridge!`; AsapOptim `SectionVariable.iglobal::Float64` and `f_axial_new` pullback early-return; AsapOptim missing `FreeFree` k method.

## Verification strategy
1. Phase 0 characterization suite is the master regression: every later phase must reproduce pinned numerics (displacements/reactions to 1e-12 where formulations are algebraically identical; documented tolerance where formulations legitimately improve, e.g. GravityLoad).
2. Parity tests: `assemble_K(cache, extract_state(model)) ≈ assemble_K!(cache)`; `solve` pure ≡ `solve!` fast.
3. AD: Zygote gradients vs FiniteDifferences on truss + frame + VariableElement problems; later DifferentiationInterface matrix (Enzyme/Mooncake) on the pure path.
4. Physics property tests: Mz′ = Vy, Vy′ = −w (finite diff along member); zero moment at released ends; semi-rigid limits reproduce release oracles; spring support reaction = −k·u.
5. Performance: `@allocated == 0` for `assemble_K!` + scalar force evaluators post-process!; benchmark suite (BenchmarkTools) comparing old vs new on a representative full-building-scale model.
