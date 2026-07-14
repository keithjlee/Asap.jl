# Asap.jl

Structural analysis (direct-stiffness FEA) for trusses and frames, plus a self-contained Force Density Method (FDM) form-finding subsystem and parametric structure generators. Pure Julia. Hard deps: LinearAlgebra, SparseArrays, StaticArrays; weak deps (extensions): ChainRulesCore, Enzyme, Mooncake.

**This is the v1.0 core** (`main` since 2026-07). The legacy v0.2.x lives on in the registry and the `v0.2.2` tag; `docs/MIGRATION_v1.md` maps the old API to the new one, `docs/MODERNIZATION.md` records the architecture rationale, `PROGRESS.md` the build-out.

## Commands

```bash
julia --project=. -e 'using Pkg; Pkg.test()'          # full suite (~2,600 tests)
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # first-time setup
julia --project=. test/readme_examples.jl             # every README code block, executed
```

Julia 1.10 is the floor; CI runs 1.10 + latest on Ubuntu/macOS.

## Architecture (three layers)

1. **Definition** — `Model{T}`: `Node`s (always 6 DOF slots), elements, loads, `NodalSpring`s. Pure user data.
2. **Analysis structure** — `AnalysisCache{T}` built once by `process!`: free/fixed/**inactive** DOF partition, element groups with plain-data mirrors, frozen sparsity + nzmap/scatter, cached factorization.
3. **Results** — `LinearResults{T}`; read via accessors: `displacement(res, node)`, `reaction(res, node)`, `element_forces(res, el)`, `axial_force(res, el)`. Results are NOT mutated onto structs.

Include order in `src/Asap.jl` is the dependency order: `Materials/` → `Nodes/` → `Elements/` (end_conditions, kernels, frame, variable) → `Springs/` → `Loads/` → `Model/` → `Analysis/` + `Results/` → `Model/utilities.jl` → `FDM/` → `Generation/` → `Geometry/` → `ShowMethods.jl` (last, so it covers every type).

Key mechanisms:
- **Unified truss/frame** — one `Model`; `TrussElement` and `FrameElement` mix freely. DOF activity (free/fixed/inactive) retires DOF blocks nothing stiffens (e.g. truss-node rotations) instead of leaving singular modes.
- **End conditions** — `EndSprings`/`EndConditions` (Monforton–Wu fixity factors). Release symbols (`:fixedfixed`, `:fixedfree`, `:freefixed`, `:freefree`, `:joist`) are exact k∈{0,∞} limits; finite k = semi-rigid, differentiable. NOTE: `:joist` keeps torsion rigid, so its rotation blocks stay ACTIVE — nodes held only by joists are mechanisms (user must restrain them).
- **`VariableElement`** — segment chains via internal DOFs (6 per interior joint), replacing the old BridgeElement shatter.
- **Loads/FEFs** — `DistributedLoad` (piecewise-linear canonical form; `LineLoad`/`TrapezoidLoad` lower to it), `PointLoad`, `PointMoment`, `SelfWeight`. One Gauss-3 FEF engine + generic condensation (`condense_fef`) replaces per-release formula catalogs.
- **Force/displacement recovery** — `internal_forces(model, el)` → `ElementForceState` with zero-allocation evaluators (`moment_z(state, t)`, …) and `InternalForces` dense sampling (stations include both sides of point actions). Displacement recovery integrates the exact curvature; start rotations come from far-end compatibility (member-side, NOT node rotations — matters at released/semi-rigid ends).
- **Cases/combos** — `solve_cases!` (multi-RHS, one factorization), `LoadCombination`, `envelope`.
- **AD-first** — pure path `extract_state`/`ModelState`/`solve` shares kernels with the in-place path; `ext/AsapChainRulesExt.jl` holds the rrules; Mooncake/Enzyme bridged via their extensions (Enzyme needs `set_runtime_activity` + `function_annotation = Const`). NEVER read mutable-struct fields inside differentiated closures — use the plain-data mirrors (see `docs/MODERNIZATION.md` and AsapOptim's `docs/`).
- **Generators** (absorbed from AsapToolkit): `Warren2D`, `Pratt2D`, `BakerTruss`, `TrussFrame`, `SpaceFrame`, `SpaceFrameBeam`, `Frame`, `GridFrame`, `GridNetwork`, ground structures + `to_truss`/`to_frame`; plot-prep `Geo`/`ElementDisplacements`; `to_network` FDM translation. Generator internals are still legacy-style Float64 (deferred polish).

## Conventions

- Parametric `{T}` core; Float64 by default. Units-agnostic — **no Unitful, ever**.
- **Documentation standard**: every function AND type gets an extensive docstring (fields, engineering meaning of symbols, dimensions in brackets like `[force/length²]` — never named units); every user-facing type gets a `Base.show(io, ::MIME"text/plain", x)` in `src/ShowMethods.jl` with aligned `symbol = value (plain-language explanation)` lines. Docs land WITH the code, never as a cleanup pass.
- Style: explicit, readable, directory-per-concept, exports in `src/Asap.jl`.

## Testing

- `test/characterization/` pins exact numerics from legacy v0.2.x and the DiffAnalysis_2024 publication env — regenerate fixtures only intentionally (a fixture diff is a behavior change).
- Deliberate deviations from legacy are documented in `docs/MODERNIZATION.md` and asserted where they occur.
- AD: Zygote/Mooncake/Enzyme vs FiniteDifferences gradient checks in-suite.

## Releasing

Bump `version` in Project.toml on `main`, then comment `@JuliaRegistrator register` on the release commit. TagBot creates the git tag + GitHub release automatically after the General registry PR merges. Never tag manually.

## Ecosystem (lockstep dependents)

- `../AsapOptim` — differentiable optimization layer (variables/params/objectives over the pure path). Registered in General. Examples + AD benchmarks in its `docs/` and `examples/`.
- `../AsapToolkit` — now only: AISC `SteelSections`, `AsapSections` polygon geometry, JSON/Grasshopper IO. Unregistered.
- `../AsapHarmonics` — connection analysis; light coupling.
- `../DemandTransport2` — optimal transport + differentiable structural analysis research; consumes AsapOptim's API; runs Julia 1.11 (pinned Makie stack).
- Tasha's fork (github.com/natashahirt/Asap.jl) — reference for TributaryLoad math; do not merge (Unitful dep).
