# Asap.jl

Structural analysis (direct-stiffness FEA) package for trusses and frames, plus a self-contained Force Density Method (FDM) form-finding subsystem. Pure Julia; deps are stdlib only (LinearAlgebra, SparseArrays).

**A major v1.0 modernization is underway — see `docs/MODERNIZATION.md` for the approved architecture and phased roadmap before making structural changes.** Until Phase 1 lands, the code described below is the *legacy* (v0.2.x) architecture; characterization tests in `test/` pin its numerics and are the regression contract for the rewrite.

## Commands

```bash
julia --project=. -e 'using Pkg; Pkg.test()'          # run test suite
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # first-time setup
```

Tests validate against textbook examples (Kassimali, *Matrix Analysis of Structures 2e*). Characterization fixtures live in `test/characterization/` — regenerate only intentionally via `julia --project=. test/characterization/generate_fixtures.jl` (they pin exact numerics; a diff in fixtures is a behavior change).

## Architecture (legacy v0.2.x)

Include order in `src/Asap.jl` is the dependency order:

- `Materials_Sections/` — `Material`, `Section` (frame: A, E, G, Ix, Iy, J, ρ), `TrussSection` (A, E, ρ)
- `Nodes/` — `Node` (6 DOF) / `TrussNode` (3 DOF); fixity via symbols (`:fixed`, `:free`, `:pinned`, `:xfree`, …)
- `Elements/` — `Element{R<:Release}`, `TrussElement`, `BridgeElement`. Releases are **type parameters** (`FixedFixed`, `FixedFree`, `FreeFixed`, `FreeFree`, `Joist`) selecting closed-form 12×12 stiffness matrices in `K.jl`. Transformation matrices in `R.jl`; `Ψ` is the roll angle.
- `Loads/` — `NodeForce`, `NodeMoment`, `LineLoad` (uniform only), `GravityLoad`, `PointLoad`; fixed-end forces in `fixedEndForces.jl` with per-release corrections
- `Model/` — `Model` / `TrussModel` (fully separate types; truss and frame **cannot mix**). Pipeline: `process!` (IDs → DOF indices → loads → sparse K assembly) then `solve!` (free-DOF sparse solve → post-process). `bridgeprocessing.jl` "shatters" host elements for BridgeElements (fragile; being replaced).
- `FDM/` — independent force-density subsystem with its own Node/Element/Load/solve!.

Truss vs frame is duplicated everywhere (separate `create_S!`, `populate_DOF_indices!`, `solve!`) — a deliberate target of the modernization.

## Conventions

- Everything is Float64 (legacy; v1.0 goes parametric).
- Units-agnostic: user is responsible for consistency (tests use kips/in, kips/ft, kN/m). **No Unitful — ever.**
- Global axes: `globalX/Y/Z` constants. Local coordinate system stored per element as `LCS`.
- Results are mutated onto structs (`node.displacement`, `element.forces` = 12-vector of local end forces; axial force = `forces[7]` frame / `forces[2]` truss).
- Style: explicit, readable, directory-per-concept, exported names in `src/Asap.jl`.

## Known legacy bugs (documented, fixed in modernization — do not silently "fix" without updating characterization tests)

- `q_local(::GravityLoad)` references nonexistent `load.element.ρ` — GravityLoad FEFs have never worked (`test/characterization/` has a `@test_broken`).
- `release!` is exported but undefined; `copy(::TrussModel)` references nonexistent `element.nodeIDs`.
- `shatterReleaseDict` in `bridgeprocessing.jl` has a duplicate `Element{FreeFixed}` key; stale second implementation `process_bridge!` is dead code.
- `scratch.jl` uses an outdated `Section` argument order.

## Ecosystem (lockstep dependents — breaking changes here must be coordinated)

- `../AsapOptim` — differentiable optimization layer (Zygote). Pinned `Asap = "0.2"`. Re-implements assembly functionally; reads `model.S` CSC internals (`all_inz`) and keeps verbatim copies of the local stiffness matrices — these must stay numerically identical to `Elements/K.jl`.
- `../AsapToolkit` — utilities; owns `InternalForces` shear/moment diagrams today (moving into core in Phase 3). Heavy direct field access; treats field names as public API. Note its `.My` is actually the moment about local z.
- `../AsapHarmonics` — connection analysis; light coupling (`connectivity`, `axial_force`, `node.reaction`, `element.LCS`).
- Tasha's fork (github.com/natashahirt/Asap.jl) — reference for load features (TributaryLoad math); do not merge (Unitful dep).
