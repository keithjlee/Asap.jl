# Asap v1.0 Modernization — Progress

> Live checklist, updated as work lands. Plan: `docs/MODERNIZATION.md`. Branch: `modernization/phase-0` (main untouched).

## Phase 0 — Hygiene & characterization ✅ COMPLETE

- [x] CLAUDE.md for Asap, AsapOptim, AsapToolkit, AsapHarmonics
- [x] Approved plan committed as `docs/MODERNIZATION.md`
- [x] CI workflow (Julia 1.10 LTS + latest, ubuntu + macos); julia compat → 1.10
- [x] Characterization fixtures pinning exact v0.2.2 numerics (94): Kassimali textbook models, portal frames × all 5 releases, per-release q()/q_local()/local_K/global_K/R oracles on a skewed element
- [x] GravityLoad bug documented as `@test_broken`
- [x] AsapToolkit InternalForces oracles (84; joist point-load combo skipped — Toolkit catalog gap)
- [x] DiffAnalysis_2024 publication golden fixtures (16): all 4 paper structures serialized + rebuilt + verified with plain Asap (0.2.1 ≡ 0.2.2 at rtol 1e-9); Zygote gradients at x_init for S4.1/S4.2 (Phase 5a AD oracle)
- [x] DemandTransport2 added to migration scope (Phase 5d) + downstream AD contract documented
- [x] Documentation standard codified (FormaSlab-style docstrings + ShowMethods.jl)
- [x] Branch policy: main = stable v0.2.x; all work on `modernization/*`

**Suite status: 6 legacy + 113 characterization passing, 1 expected broken (GravityLoad).**

### Notable findings during Phase 1 (see commit messages for detail)
- Legacy `q_local(::PointLoad)` splits axial FEF 50/50 regardless of position — physics bug, replaced by consistent lever-rule distribution (documented deviation in docs/MODERNIZATION.md).
- DOF activity must be decided per node rotation *block* in global coordinates (Λ mixes torsion with bending rotations) — encoded in `dof_signature`.
- Legacy `planarize!(model, :XY)` also zeroes element Ψ — ported.
- Legacy Asap 0.2.1 (publication) ≡ 0.2.2 forward results verified at rtol 1e-9.

## Phase 1 — Core rework ✅ COMPLETE

Order is bottom-up; each layer validated against Phase 0 oracles before the next.

- [x] `Material{T}` + `AbstractSection{T}` accessor contract (EA/EIx/EIy/GJ/ρA) + `Section{T}` + `RigiditySection{T}` — with full docstrings + show methods
- [x] `Node{T}` (SVector position, 6-DOF fixity; TrussNode absorbed)
- [x] `EndSprings{T}`/`EndConditions{T}` + release symbol map
- [x] Element kernels: `local_frame` (transformation), fixity-factor stiffness (Monforton–Wu), validated against pinned `local_K`/`R` oracles at release limits
- [x] Element structs: `FrameElement{T,S}`, `TrussElement{T,S}` + element interface (`dof_signature`, `ndofs`, …)
- [x] `NodalSpring{T}` (applicative spring supports)
- [x] DOF layer: activity accumulation, free/fixed/inactive partition, internal-DOF blocks (internal blocks land with VariableElement)
- [x] Symbolic assembly (element groups, frozen pattern, nzmap; scatter matrix lands with pure path)
- [x] In-place numeric assembly + solve pipeline (`process!`/`solve!`, cached factorization)
- [x] `LinearResults{T}` + accessors + element-wise reactions
- [x] Pure functional path (`ModelState`, `assemble_K`, `solve`) + `ext/AsapChainRulesExt.jl` (3 rrules)
- [x] Parity tests (pure ≡ in-place) + Zygote-vs-FiniteDiff gradient checks
- [x] `VariableElement{T}` via internal DOFs (bridgeprocessing deletion lands with the module flip)
- [x] Characterization suite green on new core — **module flipped**: legacy deleted, Asap IS the v1.0 core, 2347 tests green, version 1.0.0-DEV

## Phase 2 — Loads ✅ COMPLETE

- [x] Gauss/Hermite FEF engine (exact for piecewise-linear intensities)
- [x] Canonical `DistributedLoad` + `LineLoad`/`PointLoad`/`SelfWeight`; case tags on all loads
- [x] Generic end-spring FEF condensation (replaces per-release q() catalog)
- [x] SelfWeight from ρA(section) — GravityLoad bug structurally impossible
- [x] `TrapezoidLoad` convenience constructor + `PointMoment` (derivative-Hermite kernel, classic-table verified)
- [x] Trapezoid FEFs cross-checked against exact polynomial-integration oracles + classic triangle tables
## Phase 3 — Force recovery (priority) ✅ COMPLETE

- [x] `LoadTrace` (merged local loads, jump-exact doubled breakpoints, cumulative ∫w / ∫w·s)
- [x] `ElementForceState` + zero-allocation scalar evaluators (`axial_force`/`shear_y`/`shear_z`/`torsion`/`moment_y`/`moment_z` at fraction t)
- [x] `InternalForces` dense sampling (breakpoints + both sides of point actions) — axis-correct names
- [x] Displacement recovery (`local_displacements`: exact double integration, classic deflection formulas verified)
- [x] VariableElement composite dispatch (`moment_z(model, el, t)` unified queries; interior continuity)
- [x] Validated against pinned AsapToolkit oracles across all portal releases under the documented rename map (tk.My→Mz, tk.Mz→−My, tk.P→N)
- [x] Property tests: dMz/dx=Vy, dMy/dx=−Vz, dVy/dx=w, exact jumps at point forces/moments
## Phase 4 — Cases/combos/envelopes ✅ COMPLETE

- [x] `solve_cases!`: per-case solves against ONE assembly + ONE factorization (multi-RHS)
- [x] `LoadCombination` + `combine`: factored results by pure superposition — verified ≡ brute-force factored solves (u, reactions, forces, compliance, recovery diagrams)
- [x] Combination-aware recovery (`internal_forces(...; results, factors)`)
- [x] `envelope`: station-wise extrema over combinations at common stations (all-cases-on station set so no combo's jumps are missed) — verified ≡ brute-force extrema
- [x] Single-factorization reuse proven in tests
## Phase 5 — Ecosystem migration 🔄

### 5a AsapOptim ✅ COMPLETE (branch `asap-v1` in ../AsapOptim)
- [x] Rewritten as a ~600-line Zygote-free design-vector layer on the core pure path (was 2400 lines of parallel implementation + hand-written adjoints)
- [x] Unified `OptParams` (Truss/Frame aliases); Spatial/Area/Coupled variables with sparse scatter compilation
- [x] `solve_structure`/`axial_force`/`axial_stress`/`compliance`/`GeometricProperties`/`updatemodel`
- [x] First-ever AsapOptim test suite: 28 tests, gradcheck vs FiniteDifferences, design-eval ≡ direct solve
- [x] Legacy preserved unloaded in `legacy_v0/`
- [x] FDM network optimization path (QVariable/NetworkOptParams/solve_network — forward parity 1e-10, gradcheck green; 32 AsapOptim tests)
- [ ] Deferred: SectionVariable (RigiditySection parameterization is the v1.0 idiom)

### 5b AsapToolkit ✅ PORTED (branch `asap-v1`) — visual verification done (2026-07-14)
- [x] ForceAnalysis/ deleted (absorbed into core); generators/Geometry/IO/SteelSections/FDM-translations ported
- [x] Frame generator: BridgeElement joists → explicit interior nodes; releases now actually applied (legacy passed them as the id positional — silently ignored!)
- [x] ElementDisplacements reimplemented on core displacement recovery
- [x] First test suite: 58 tests (generator smoke + equilibrium, displaced shape vs core, section bridges)
- [x] Side-by-side plot verification vs legacy stack (Warren2D, beam diagrams, Frame generator, SpaceFrame): all numerics identical to 7+ significant digits; comparison page published as a Claude artifact ("Asap v1.0 — visual verification vs legacy"). Frame case confirms the release-as-id legacy bug numerically: v1.0 forced to all-fixedfixed reproduces legacy to 8+ digits; with the requested releases, peak displacement differs ~17%. Also confirmed legacy Toolkit main's `InternalForces(element, model)` and `ElementDisplacements` are dead code vs Asap 0.2.2 (read `element.release` / deleted `release2DOF`/`DLineLoad`/`DPointLoad` tables); only the `etype2DOF`-keyed methods work.
- [x] Keith eyeballed the comparison plots (2026-07-14): approved, 5b closed.

### AsapToolkit → Asap absorption (2026-07-14, Keith-approved clean break) ✅
- [x] Moved into Asap core (zero new dependencies): `Generation/` (all generators + ground structures; Statistics `mean` inlined; Interpolations requirement removed — variable-depth `SpaceFrame`/`CornerDome` now take any callable `surface(u, v)` on [0,1]²), `Geometry/` (`Geo`/`ModelGeo`/`TrussGeo`/`NetworkGeo` + `ElementDisplacements`/`displacements`), `to_network` (FDM translation), `clear_supports!`/`element_connectivity` (`Model/utilities.jl`)
- [x] Documentation standard: struct docstrings for all moved types; `Base.show` methods (generic generator/ground-structure display + Geo/ElementDisplacements) in `ShowMethods.jl` (now included last)
- [x] Tests moved with the code: `test/newcore/test_generation.jl` (Asap suite now 2635) + generators/geo example block in README (tested via `test/readme_examples.jl`)
- [x] AsapToolkit slimmed to `SteelSections` + `AsapSections` + IO (`topologize`, `GHsave`); deps down to Asap/JSON/LinearAlgebra/Reexport/XLSX; no re-exports (clean break); suite 14 tests green; Keith's AsapSections WIP untouched
- [x] DemandTransport2 dropped its AsapToolkit dependency entirely (`TrussGeo` now from Asap); 108 tests green

### 5c AsapHarmonics ✅ PORTED (branch `asap-v1`)
- [x] Rename pass complete; first test suite (3D spherical + planar 2D signatures, finite feature vectors)

### 5d DemandTransport2 ✅ PORTED (branch `asap-v1`; package src) — examples pending
- [x] Package src ported (generators, ElementDemand, objectives, ResponsiveGeo); first test suite: 108 tests incl. gradient checks (Julia 1.11 — pinned Makie stack)
- [x] Latent bugs fixed: restore! @lift parsing; Asap/Nonconvex Model export ambiguity
- [ ] Example scripts (highrise_to_bridge, ot_example) — port interactively with Keith
- [ ] spaceframe_to_vault — network path now available; script port is interactive work

## AD verification & benchmarks (2026-07-13)
- [x] Cross-stack verification on the publication S4.2 spaceframe: new-stack gradients ≡ legacy publication stack (compliance 3.6e-13, volume bit-identical); displacements ≡ publication at 1.1e-13
- [x] Performance: solve 1.2× faster, assembly 1.5× faster, differentiable forward 2× faster; ∇compliance 2.8 ms (within 2× of the fully hand-differentiated legacy layer, fully generic); full report + reproducible scripts in `AsapOptim/docs/`
- [x] Differentiable-path optimization: batched truss assembly, plain-data mirrors, analytic truss_stiffness rule (6631 ms → 2.8 ms)
- [x] README rewritten in depth with a TESTED example per feature (`test/readme_examples.jl`)

## AD backend liberation (2026-07-13) — escape from Zygote ✅
- [x] Mooncake works (AsapMooncakeExt: @from_rrule bridge + sparse-cotangent adapter + CHOLMOD Factor tangent-free) — Zygote-accuracy, 3× less memory at scale
- [x] Enzyme works (AsapEnzymeExt: @import_rrule in __init__; needs runtime-activity + Const function annotation) — **fastest backend, beats the legacy hand-differentiated stack** (1.04 vs 1.41 ms, 512 vars)
- [x] Core rules representation-flexible; solve_free cotangent projected onto the frozen pattern (dense outer product eliminated for all backends)
- [x] Report: `AsapOptim/docs/AD_BACKENDS.md` (usage snippets, gotchas, reproduce script)

## Examples & JointVariable (2026-07-13)
- [x] NEW core capability: `ModelState.ends` override + `JointVariable` — connection stiffness as a differentiable design variable (FEF-consistency guard for element loads on jointed members)
- [x] AsapOptim examples ported + expanded, ALL executed before commit: truss-optimization1/2 (v1.0 port), vierendeel-geometry (frame-action C·V optimization), joint-stiffness (precast frame, pinned bases — drift lands exactly on H/400 at 1.1% of the uniform stiffness budget), ad_backends/ (Zygote/Mooncake/Enzyme, one file each, identical gradients)

## :joist soundness audit (2026-07-14, Keith-requested) ✅
- [x] Verified k(:joist) == axial ⊕ torsion exactly (rank 2; Monforton–Wu k→0 limit to 2.7e-15; k→∞ → fixedfixed to 8.6e-13)
- [x] Statics: UDL joist delivers wL/2 end shears, zero end/reaction moments; internal recovery gives wL²/8 midspan moment, zero end moments; torsion path exact (TL/GJ)
- [x] FOUND + FIXED: displaced-shape recovery used NODE rotations for u0 — wrong at released/semi-rigid START ends (joist midspan recovered wL⁴/128EI instead of 5wL⁴/384EI). Fixed via far-end-compatibility start rotations; applies to frame/truss/VariableElement states; +6 regression tests (suite 2641)
- [x] Documented modeling caveat: :joist keeps the whole rotation block ACTIVE (torsion rigid ⇒ block can't be retired since Λ mixes torsion+bending in global DOFs) — bending rotations and transverse translations at joist-ONLY nodes are mechanisms; solver raises ZeroPivotException/PosDefException. User must restrain rotations at nodes held only by joists.

## Deferred / follow-ups
- SectionVariable → RigiditySection parameterization idiom
- FDM subsystem parametrization in core
- Parametrize the generators (`{T}` structs, concrete fields) — internals are still legacy-style Float64; API-ported and verified, so this is polish, not correctness
- Geometric nonlinearity (`geometric_stiffness` hook ready), dynamics (`mass` hook ready)
- ✅ ~~Enzyme/Mooncake/DifferentiationInterface matrix on the pure path~~ (done 2026-07-13 — see AD backend liberation above)
- Batch _element_lengths / frame-group assembly for further gradient speed (notes in AsapOptim/docs/AD_VERIFICATION_AND_BENCHMARKS.md)
- Docs build (Documenter.jl) from the extensive docstrings
