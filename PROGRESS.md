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

## Phase 2 — Loads 🔄 (largely landed during Phase 1)

- [x] Gauss/Hermite FEF engine (exact for piecewise-linear intensities)
- [x] Canonical `DistributedLoad` + `LineLoad`/`PointLoad`/`SelfWeight`; case tags on all loads
- [x] Generic end-spring FEF condensation (replaces per-release q() catalog)
- [x] SelfWeight from ρA(section) — GravityLoad bug structurally impossible
- [ ] `TrapezoidLoad` convenience constructor + `PointMoment`
- [ ] Trapezoid FEFs cross-checked against closed-form oracles
## Phase 3 — Force recovery (priority) ⬜
## Phase 4 — Cases/combos/envelopes ⬜
## Phase 5 — Ecosystem migration (5a Optim / 5b Toolkit / 5c Harmonics / 5d DemandTransport2) ⬜
