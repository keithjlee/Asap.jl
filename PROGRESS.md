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

## Phase 1 — Core rework 🔄 IN PROGRESS

Order is bottom-up; each layer validated against Phase 0 oracles before the next.

- [x] `Material{T}` + `AbstractSection{T}` accessor contract (EA/EIx/EIy/GJ/ρA) + `Section{T}` + `RigiditySection{T}` — with full docstrings + show methods
- [x] `Node{T}` (SVector position, 6-DOF fixity; TrussNode absorbed)
- [x] `EndSprings{T}`/`EndConditions{T}` + release symbol map
- [ ] Element kernels: `local_frame` (transformation), fixity-factor stiffness (Monforton–Wu), validated against pinned `local_K`/`R` oracles at release limits
- [ ] Element structs: `FrameElement{T,S}`, `TrussElement{T,S}` + element interface (`dof_signature`, `ndofs`, …)
- [ ] `NodalSpring{T}` (applicative spring supports)
- [ ] DOF layer: activity accumulation, free/fixed/inactive partition, internal-DOF blocks
- [ ] Symbolic assembly (element groups, frozen pattern, nzmap + scatter)
- [ ] In-place numeric assembly + solve pipeline (`process!`/`solve!`, cached factorization)
- [ ] `LinearResults{T}` + accessors + element-wise reactions
- [ ] Pure functional path (`ModelState`, `assemble_K`, `solve`) + `ext/AsapChainRulesExt.jl` (3 rrules)
- [ ] Parity tests (pure ≡ in-place) + Zygote-vs-FiniteDiff gradient checks
- [ ] `VariableElement{T,S}` via internal DOFs; delete bridgeprocessing
- [ ] Characterization suite green on new core

## Phase 2 — Loads ⬜
## Phase 3 — Force recovery (priority) ⬜
## Phase 4 — Cases/combos/envelopes ⬜
## Phase 5 — Ecosystem migration (5a Optim / 5b Toolkit / 5c Harmonics / 5d DemandTransport2) ⬜
