# Shell elements in Asap — implementation plan

> **Status: APPROVED 2026-11-23 (Keith). In progress.** Branch `shells` (off `main` @ `1059e4f`).
> Target: Asap **1.3.0** (additive; no public-API break). Lockstep suites run at every phase boundary.
> Companion files: `DERIVATION.md` (S0 dossier, formulations and benchmark values) and the
> progress tracker below. Keep that tracker current; it is the resume point.

## Resume protocol (read first when picking this back up)

1. `git switch shells && git log --oneline main..` shows what has landed. Each phase ends in a
   commit whose subject starts with `shells(Sn):`.
2. Read the **Progress tracker** below. The first unchecked item is the next task. The
   "Handoff notes" under each phase record state that the code doesn't show (half-done work,
   decisions, gotchas).
3. Every phase is **independently shippable to the branch**: the full Asap suite stays green at
   each phase commit. Stopping mid-phase means leaving a `WIP` note in the handoff, never a
   red suite on the branch.
4. Baseline numbers: Asap suite 2 641+ tests, about 4–6 min on 1.13 (see the S1 handoff for the
   exact count at branch-off).

## Progress tracker

| Phase | Status | Commit |
|---|---|---|
| S0 derivation dossier | ✅ membrane transcribed from ITW report; MITC4 written, primary-PDF recheck open | — |
| S1 n-node analysis-layer generalization | ✅ bitwise parity with main; lockstep green | `e5a8ce2` |
| S2 `ShellSection` + ITW drilling membrane | ✅ enhanced variant (DERIVATION §2.1) — Keith OK'd 2026-11-23 | `4e9677f` |
| S3 MITC4 plate → full `ShellElement` | ✅ plates + Scordelis–Lo converge (DERIVATION §3) | `2b67b1e` |
| S4 area loads, resultants, `section_cut` | ✅ | `da60f78` |
| S5a wall mesher + edge loads + 2-node link element | ✅ | `bc682a4` |
| S5b rigid diaphragm constraints | ✅ exact transformation; torsion hand-calc to 1e-9 | `deddfa7` |
| S6 pure path / AD + AsapOptim thickness variable | ✅ design-dependent shell loads; `AsapShellADExt`; AsapOptim `ThicknessVariable` (branch `shells`) | (see `git log`) |
| S6b feature parity (Keith, 2026-11-23 — before release) | ✅ P1–P7; lockstep green (Asap 3067), docs + README pass | `0826a8e`, `e7e85da`; AsapOptim `0dc0ae8` |
| S6c examples for ALL scenarios + Enzyme across Asap (Keith, 2026-11-23 evening) | ✅ walkthroughs + Typst reports + Enzyme matrix (1.13 + 1.12, 18/18) | see below |
| S6d reverse-mode performance (2026-11-23/24) | ✅ AsapFrameADExt; gather refactor REVERTED (broke Enzyme); linear-cost RuleConfig group rules + split shell rules (vs v1.2 frames 8.8–11.7×; vs kernel rules only frames 2.4–6.9×, shells 2.4–4.7×); bitwise primal vs v1.2; Enzyme 36/36 | release notes "Performance" |
| Release 1.3.0 (docs pages, release notes, lockstep, registrator) | ⏳ everything ready; **awaiting Keith: merge + bump + push + register** | |
| S7 triangles, mesh import, warping correction (post-1.3) | ☐ | |

### S6c — walkthroughs for every scenario + Enzyme verification (Keith: before release)

Keith: "full walkthrough example docs for all structural scenarios in the examples folder:
frames, trusses, 2D/3D, varying joint stiffnesses, spring supports, shape vs size optimization,
etc." and "verify that Enzyme can now work across Asap".

Layout: `Asap/examples/` is ONE environment (Project.toml moved up from `shells/`).
`figures.jl` is shared (line-element plots added; `set_figdir(@__DIR__)` per folder).
Folders: `analysis/`, `shells/`, `optimization/`, `ad/`.
- [x] E1 restructure: shared env and figures; shells walkthroughs re-pointed and re-run
- [x] E2 `analysis/01_truss_2d` — planar truss by hand and by generator: supports, loads,
      axial forces, reactions, deflected shape
- [x] E3 `analysis/02_spaceframe_3d` — 3D space truss: self-weight, reactions, force distribution
- [x] E4 `analysis/03_frame_2d` — portal frame: line/point loads, releases, internal-force
      diagrams, deflected shape
- [x] E5 `analysis/04_frame_3d_building` — 3D building frame: load cases, combinations,
      envelopes, drift
- [x] E6 `analysis/05_joints_springs_variable` — semi-rigid connections (EndSprings sweep),
      spring supports (NodalSpring, soil), VariableElement (haunched girder)
- [x] E7 `optimization/01_truss_size` — AreaVariables: compliance with a volume limit
      (optimality criteria), gradients checked
- [x] E8 `optimization/02_truss_shape` — SpatialVariables (coupled/symmetric): shape vs size,
      what each buys
- [x] E9 `optimization/03_frame_joints` — JointVariable: cheapest connection stiffness
      meeting drift
- [x] E10 `ad/enzyme_matrix.jl` — Enzyme (reverse and forward) vs Zygote/FD on: truss area,
      truss shape, large spaceframe, frame with element loads, semi-rigid joints, springs,
      VariableElement, shells (thickness), diaphragm. One process per case under `timeout`;
      Julia 1.13 (+1.11 if time); results table in `ad/ENZYME_STATUS.md` and the memory note
- [x] E11 index README, docs guide page "Examples", lockstep, commit

### S6b — feature parity checklist (Keith: "full feature parity for shells before 1.3")

Audit 2026-11-23 of user-facing paths on models containing shells:
- [x] P1 `solution_tangents` (AsapOptim implicit diff) with shells: ThicknessVariable,
      SpatialVariables on shell nodes, design-dependent shell loads, rigid diaphragms
- [x] P2 plot extraction: `Geo(model)`/`ModelGeo` currently break on mixed models (assume line
      elements); line part must skip shells, plus a new `ShellGeo` (quads, displaced
      positions, per-element and nodal-averaged resultant fields for contouring)
- [x] P3 `shell_envelope` — resultant extrema over load combinations (frames have `envelope`)
- [x] P4 `element_connectivity` generic over N-node elements (currently via line-only
      `connectivity`)
- [x] P5 `Model` show counts shells/links/diaphragms
- [x] P6 AD backends: Mooncake through shells (bridge the shell rrules if needed); Enzyme
      status documented
- [x] P7 walkthroughs `examples/shells/` (public; own env: Asap path + CairoMakie + kjlMakie):
      shift+enter friendly, top-to-bottom runnable — shells standalone (wall, slab), mixed
      frame+wall building with diaphragm, precast joints, thickness gradients/optimization;
      `figures.jl` figure generators (mesh, deformed shape, resultant contours, pier forces)
      P6: Mooncake needed `AsapShellMooncakeExt` (its own rule derivation overflows on the
      shell kernels) plus cotangent adapters (nested `Tangent` in, static-array rdata out);
      Mooncake ≡ Zygote to 1e-16. CORRECTION (2026-11-23 evening): an earlier "Mooncake faster,
      13 vs 49 ms" compared a 12- with a 29-element model. On the same models Mooncake is the
      SLOWEST (96 el: Zygote 51, Enzyme 32, Mooncake 84 ms). Enzyme REVERSE works via
      `AsapShellEnzymeExt` (`@import_rrule` of the shell rules; native Enzyme fails in type
      analysis on the large SMatrix registers): ≡ Zygote to 1e-16 and the FASTEST.
      Enzyme forward mode is not bridged for shells; use ForwardDiff. Benchmark script:
      `examples/ad/ad_backends_shells.jl`.
- [x] P9 Bathe & Dvorkin 1985 primary-source check (`test_shell_mitc4_bd85.jl`;
      DERIVATION §3).

### Handoff notes

- **S2 (done):** `ShellElement{T,4,S}`, `ShellSection`/`ShellRigiditySection` (+ modifier
  constructor), kernels in `src/Elements/kernels/shell.jl`. The element is currently
  **membrane-only**: `_plate_stiffness` returns zeros, so out-of-plane DOFs must be
  restrained (tests planarize) until S3 lands MITC4 there. Kernel performance: about 5 µs and 0
  allocations per 24×24 global stiffness; first-call compile about 8 s.
  - **Gotchas:** (1) `SVector{3}(::SVector)` has no working Zygote pullback; use `_sv3`. (2)
    Closures capturing reassigned locals get boxed (294 allocations); keep one function barrier
    per quadrature term. (3) `kron` and 24×24 0/1 embedding products are slow and Zygote-hostile;
    use the static-index gather (`_shell_local`) and the reshape-based blockwise rotation
    (`_rotate_shell`). (4) With a drilling element, "clamped" must fix ψ too.
  - ForwardDiff through the kernel works but its first compile is about 50 s (Dual SMatrix
    24×24); revisit in S6.
  - Benchmark oracles live in `test/newcore/shell_test_utils.jl` plus `test_shell_membrane.jl`.
- **S4 (done):** `AreaLoad` (global/local, pressure shorthand) and `SelfWeight` on shells
  lower to CONSISTENT loads, including drilling moments from in-plane traction and the condensed
  bubble. `shell_resultants(model, el; at, results, factors)` gives N/M/Q in local axes, and
  `section_cut(model, elements, cut_nodes; about, axes, results)` is built from nodal forces,
  so it is exact for statics. Shell loads are FROZEN in the pure path (static, like
  `VariableElement` loads) until S6.
  - **Gotcha:** a uniform LINE traction on a drilling-element edge is consistent only with
    end moments ∓p·l²/12 per edge (−at J, +at K). Nodal forces alone are not. S5a should ship
    an edge-load helper that adds them.
- **S5a (done):** `mesh_wall` (structured, with openings and hard lines; `reuse` connects
  walls to walls and frames), `wall_cut` (pier/spandrel free bodies for `section_cut`),
  `edge_loads` (consistent line traction including the drilling end moments), and
  `LinkElement` (zero-length 6-spring, enforced at `process!`). Verified: a beam framing
  in-plane into a stiff wall reaches PL³/3EI within 2% (drilling moment transfer); a stiff
  link equals the monolithic wall to 1e-6; pier shears sum exactly to the storey shear.
- **S5b (done):** `RigidDiaphragm(nodes; master, axis, id)` plus the `Model(...; diaphragms)`
  keyword (auto-masters join `model.nodes`). `ConstraintMap` lives in the cache and every solve
  path has a constrained branch: `solve!`, `solve_cases!`, pure `solve`/`compliance`.
  Unconstrained branches are untouched, and existing models are still bitwise identical to
  `main`. In the pure path C (the offsets) is frozen at the last `solve!`.
  - **Note:** `Model` gained a field (`diaphragms`, before `cache`). Downstream code reads
    fields by name, so nothing broke, but AsapOptim should be aware.
- **S6 (done):**
  - `FEFTable` shell-load fields; `assemble_Pf` re-evaluates `AreaLoad`/`SelfWeight` on shells
    per state.
  - Zygote ≡ finite differences for positions and thickness. A `:frozen` thickness gradient
    differs from the design-dependent one, because self-weight follows t.
  - **Perf:** Zygote tracing through the 24×24 kernels was 796 ms per gradient (29 elements;
    primal 1.9 ms) with a 128 s compile. `ext/AsapShellADExt.jl` adds `rrule`s on
    `_shell_K_plain`/`_shell_fef_plain` that compute the local Jacobian with ForwardDiff
    (chunk 4: 85 s compile, 49 ms gradient; chunk 12: 155 s / 52 ms; chunk 2: 75 s / 101 ms).
    Remaining cost: first-call compile.
  - **Follow-ups:**
    - batch shell kernels across elements (cf. the deferred batched-frame work);
    - an analytic membrane/plate pullback if compile time matters;
    - make `solution_tangents` shell-aware (it currently refuses thickness variables and
      diaphragms).
  - AsapOptim (branch `shells`):
    - `ThicknessVariable(elements, value, lb, ub)` covers whole walls;
    - `shell_volume(x, p)`;
    - generic `i1`/`i2`;
    - **two bugs fixed in passing:** `solve_structure` would have frozen shell loads (new
      `_design_F` predicate), and it bypassed rigid diaphragms (now through
      `Asap._constrained_solve`).
- **Release (pending Keith):**
  - docs pages `guide/shells.md` and `api/shells.md`, `docs/releases/v1.3.0.md`, the README
    section, and `test/readme_examples.jl` (its stale pre-ecosystem path is fixed);
  - steps: merge `shells` to `main` in Asap and AsapOptim, bump to 1.3.0 (AsapOptim compat
    `Asap = "1.3"` if it needs the new API, plus its own version), push, docs CI, and
    `@JuliaRegistrator register`.
- **S3 (old note, done):** implement `_plate_stiffness(xy, section)` as MITC4 (DERIVATION §3) in the
  12×12 (w, θx, θy)×4 layout, which `_shell_local` already embeds. First pull Bathe & Dvorkin
  1985 (Undermind) to confirm the tying-point formulas.

- **S0:** the ITW report PDF (UCB/SEMM-89/17, scanned, no text layer; read as page images) was
  fetched from escholarship (`qt9v61q2nj`). Bathe & Dvorkin 1985 and MacNeal & Harder 1985 are
  not yet pulled; Undermind is the route. Their numbers in `DERIVATION.md` §3 are from memory
  and are flagged.

## Decisions (Keith, 2026-11-23)

1. Plate part: **MITC4**.
2. Membrane: **ITW drilling quad** (D-type, condensed bubble).
3. Type design: **one `ShellElement{T,N,S}`**, parametric in node count N (4 now, 3 in S7); the
   kernel is selected by N. If a second formulation for the same N is ever wanted, add a
   formulation type parameter with a default, not a new type.
4. 1.3 scope includes the **mesher, link element and rigid diaphragms**.
5. External cross-check: no OpenSees/SAP today. Keith has **Karamba3D**, and OpenSees can be
   installed at validation time. Primary-literature benchmark values remain the in-suite oracle.
6. **Pure-path AD support ships in 1.3.**

## 0. Why, and what "done" means

Asap today has only line elements (`TrussElement`, `FrameElement`, `VariableElement`). Walls are
modelled as equivalent frames (the wide-column analogy). That works for slender, solid walls. It
breaks down for:

- squat walls and deep wall-beams, where shear deformation and in-plane stress fields dominate;
- walls with openings (piers and spandrels), where the stiffness and force distribution depend on
  where the openings are;
- load sharing among walls, cores and frames, which depends on relative in-plane stiffness;
- semi-rigid floor diaphragms;
- eventually, any generic curved or folded shell geometry.

**Done for 1.3** means four things:

1. a flat 4-node shell element (membrane with drilling DOF + plate bending) mixes freely with
   frames and trusses in one `Model`;
2. it passes the standard benchmark set (§5);
3. it reports stress resultants and **section-cut forces** (pier and spandrel P, V, M), the
   quantities wall design actually consumes;
4. a mesher turns a rectangular wall with rectangular openings into a conforming mesh.

Generic shells (triangles, imported meshes, warped geometry) come after, on the same machinery.

## 1. Element choice

### Options considered

| Family | Examples | Verdict |
|---|---|---|
| Equivalent frame / wide-column | today's approach | Keep as the fast option for slender walls. Not a replacement for shells. |
| Plane-stress membrane only | CST, bilinear Q4, QM6 (incompatible modes), Allman triangle, drilling quads | Enough for in-plane wall action **alone**. Plain Q4 and CST are **not acceptable**: parasitic shear makes a wall meshed one or two elements wide far too stiff in in-plane bending. |
| Flat shell = membrane ⊕ plate | MITC4 (Dvorkin–Bathe), DKQ/DKT (Batoz), DKGQ (GQ12 + DKQ, OpenSees), SAP/ETABS thin/thick shell | **Chosen.** Walls see out-of-plane loads too (wind, earth, eccentric bearing), and the same element covers slabs and diaphragms, and later generic shells. |
| Curved / degenerated-solid shells, solid bricks | MITC9, solid-shell | Overkill for building work. Revisit only for real curved-shell research. |
| Nonlinear wall macro-elements | MVLEM, SFI-MVLEM | Belongs to the deferred nonlinear phase, not here. |

### Recommendation: 4-node flat shell, assembled from two well-understood parts

- **Membrane: a quadrilateral with a true drilling rotation.** This is Ibrahimbegović–Taylor–Wilson
  (1990), built on the Hughes–Brezzi (1989) variational drilling formulation.
  - Why drilling matters *in Asap specifically*: DOF activity is decided per node **rotation
    block in global coordinates** (`dof_signature`), because Λ mixes components. A membrane
    without drilling stiffness leaves an unstiffened rotation about the element normal. For any
    wall not aligned with a global plane, that rotation cannot be retired, so it surfaces as a
    singular mode, unless we add a fictitious spring (the SAP/ETABS hack).
  - A real drilling DOF also lets a frame framing **in-plane** into a wall node transfer moment,
    which matters for coupling beams and spandrels. That is the classic "beam pinned into shell"
    modelling error.
  - Candidate refinement: add incompatible/enhanced in-plane bending modes on top of the drilling
    quad. Decided by the S2 coarse-mesh benchmark, not up front.
- **Plate bending: MITC4** (Dvorkin & Bathe 1984; Bathe & Dvorkin 1985). This is Reissner–Mindlin
  with assumed transverse-shear strains: no shear locking, valid for thin and thick plates. It is
  the most-validated quad in practice (ADINA, OpenSees `ShellMITC4`).
  - Alternative: **DKQ** (Batoz & Tahar 1982), which is Kirchhoff/thin only. It is cheaper and
    very accurate for thin slabs and walls, but gives no transverse-shear resultants and can't do
    thick plates. OpenSees `ShellDKGQ` (GQ12 membrane + DKQ) is the reference for this pairing.
  - Recommendation: MITC4, because one element is then valid across all thicknesses. Keep a
    `DKQ` plate kernel as a possible S3 spike if MITC4 under-performs on coarse slab meshes.
- **Later (S7): a 3-node companion** for mesh transitions and generic surfaces: an Allman-type
  drilling triangle ⊕ DKT or MITC3(+). Also MITC4+ (Ko, Lee & Bathe 2017) if distorted-mesh
  accuracy becomes the bottleneck.

## 2. Architecture impact (the real work)

The element *interface* is already n-node in principle: `ndofs = 6·length(nodes(el))`, and
`dof_signature` is a tuple. The **analysis layer, however, hard-codes two nodes** in these places:

| Place | Two-node assumption |
|---|---|
| `element_global_dofs(el)` | `nodeStart`/`nodeEnd`, a 12-vector |
| `stiffness(el, x1, x2)` and `_assemble_group!` | two position arguments |
| `ElementGroup` | `i1`, `i2`, `rollangles`, `endss`, `Cinc` (an nel × nnodes incidence with ±1) |
| `AnalysisCache.q_local` | `zeros(T, 12)` per segment |
| `FEFTable` / `assemble_Pf` | 12-vector per load, Λᵀ rotation |
| `ModelState` / `_group_entries` | `_position(X, i1)`, `_position(X, i2)` |
| Results | `element_forces`, `internal_forces`, `connectivity`, `Geo` |

### Plan: generalize internally, keep every public 2-node API and every field downstream reads

- Element positions come through one generic accessor: `element_positions(el, X) -> NTuple{N,SVector{3}}`.
  Line elements keep `stiffness(el, x1, x2)` as a method.
- Add a node-index matrix `inodes::Matrix{Int}` (N × nel) to `ElementGroup`, and use it for
  `element_global_dofs`. **Keep `i1`/`i2`/`Cinc` for line groups.** AsapOptim
  (`parameters.jl`, `results.jl`, `implicit.jl`) and AsapHarmonicsOptimExt read them. For shell
  groups they are empty or unused.
- Per-element "local force" storage becomes element-type-specific. Line elements keep `q_local`.
  Shells store their nodal-displacement-derived resultants lazily in results (§4), not in the cache.
- `connectivity(model)` and `Geo` stay line-element only. Shells get their own
  accessors (`shell_elements(model)`, `ShellGeo`).
- **Gate: S1 is a pure refactor.** The full Asap suite must stay bitwise-identical (characterization
  fixtures included), AsapOptim/AsapHarmonics suites stay green, and the AsapOptim
  gradient benchmarks stay unregressed.

### Performance and AD design of the kernel

- A 24×24 `SMatrix` is past the comfortable StaticArrays size (compile time, register pressure).
  Build in the element's local plane as **separate 12×12 blocks**: membrane (u, v, θz) × 4 and
  plate (w, θx, θy) × 4. Interleave them into node-wise 6×6 blocks, then rotate blockwise with one
  3×3 local frame, as frames already do with Λ. Spike this in S2: compile time, allocations, and
  `@btime` versus a frame element.
- The kernel is a pure function of (node positions, `ShellSection`), with Gauss loops over
  `SVector`s, so it is AD-transparent from day one. Per-element closures in the pure path at first.
  Batch them later together with the deferred batched-frame-kernel work (same lever, same notes).
- The local frame of a flat quad: normal = cross product of the diagonals; local x from the
  projected first edge, or from a user material-axis vector (`ShellSection`/element field, needed
  for orthotropic sections and for reporting resultants in "wall axes"). Warped quads get projected
  onto the mean plane. Warping correction is deferred to S7, and a warping-ratio warning is raised
  in S2.

## 3. Sections and materials

New `ShellSection{T} <: AbstractShellSection{T}`, mirroring the frame-section **rigidity accessor
contract**:

| Accessor | Meaning | Dimension |
|---|---|---|
| `Dm(sec)` | membrane rigidity (3×3) | [force/length] |
| `Db(sec)` | bending rigidity (3×3) | [force·length] |
| `Ds(sec)` | transverse shear rigidity (2×2, incl. κ = 5/6) | [force/length] |
| `ρt(sec)` | mass per unit area | [mass/length²] |

- `ShellSection(material, t)`: isotropic, homogeneous.
- `ShellRigiditySection(Dm, Db, Ds, ρt)`: effective or orthotropic stiffness, the analogue of
  `RigiditySection`. This covers ribbed or hollow-core decks as equivalent plates.
- **Stiffness modifiers** (f11, f22, f12, m11, m22, m12, v13, v23) as a convenience constructor, so
  code-mandated cracked stiffness (e.g. ACI 318-19 Table 6.6.3.1.1(a): walls 0.70/0.35 Ig, slabs
  0.25 Ig) is one keyword. The code-specific values live in user or kernel code, not in Asap.
- Layered/fiber sections stay out of scope; they belong to the nonlinear phase.

## 4. Loads and results

**Loads** (consistent nodal loads via Gauss quadrature, in the same spirit as the FEF engine):
- `AreaLoad(el, q; direction = :normal | global vector)`: uniform pressure, later bilinear;
- `SelfWeight` extended to shells (ρt · area);
- `EdgeLoad` (line load along an element edge) — maybe; often a frame along the edge is cleaner.
- Pure path: area loads re-evaluated per state (the area changes with node positions), riding on a
  generalized `FEFTable` (per-load N·6 vectors).

**Results** (computed from the solved displacements, never mutated onto structs):
- `shell_resultants(res, el; at = :gauss | :nodes | :center)`: Nxx, Nyy, Nxy [force/length],
  Mxx, Myy, Mxy [force], Qx, Qy [force/length], in element or material axes;
- nodal averaging and smoothing as an explicit, opt-in step (never silently);
- **`section_cut(res, elements, plane)`**: integrated P, V2, V3, T, M2, M3 over a cut through a set
  of shell (and frame) elements, about a chosen reference point. It is computed from **element
  nodal forces** (K_e u_e) on one side of the cut, which is exact for equilibrium, not from
  smoothed stresses. This is the pier/spandrel force a wall design check consumes.
- A `Pier`/`Spandrel` labelling helper on top of `section_cut`.
- Cases, combinations and envelopes: shell results superpose linearly, so `combine`/`envelope`
  extend directly.
- Show methods and docstrings land with each type (house standard).

## 5. Verification set (primary-literature reference values)

All in-suite, with reference values cited in the test file.

- **Patch tests**: membrane constant-stress and plate constant-curvature patch tests on distorted
  meshes; rigid-body modes (6 zero eigenvalues per free element); drilling-mode check.
- **Membrane**: Cook's membrane; MacNeal–Harder straight cantilever (rectangular, trapezoidal and
  parallelogram meshes; in-plane and out-of-plane tip loads); a **deep cantilever wall meshed one
  and two elements wide versus Timoshenko beam theory with shear deformation**. That last test is
  the key wall test.
- **Plate**: simply supported and clamped square plates under uniform load (Timoshenko &
  Woinowsky-Krieger, e.g. SS centre deflection 0.00406 qa⁴/D); thin-to-thick sweep (locking
  check); Morley 30° skew plate.
- **Shell**: Scordelis–Lo roof, pinched cylinder, twisted beam (warped elements), pinched
  hemisphere (MacNeal & Harder 1985 values).
- **Mixed models**: a frame beam framing in-plane into a wall (drilling transfer); a coupled wall
  with spandrels versus a published or closed-form coupled-wall solution; a wall-frame building
  where Σ reactions = Σ applied load.
- **External cross-check** (optional): the same models in OpenSees (`ShellMITC4`/`ShellDKGQ`) and,
  if available, SAP2000/ETABS.
- **AD**: gradients of compliance and a section-cut force with respect to thickness and node
  coordinates, checked against ForwardDiff, Zygote and finite differences.

## 6. Phases

| Phase | Content | Gate | Rough size |
|---|---|---|---|
| **S0** | Derivation dossier: pull the primary papers (Ibrahimbegović 1990, Hughes–Brezzi 1989, Dvorkin–Bathe 1984/85, Batoz–Tahar 1982, MacNeal–Harder 1985), write the formulations in Asap notation (local frame, B-matrices, assumed-strain tying points, drilling penalty γ). | Keith reviews the derivation | 1 day |
| **S1** | n-node generalization of the analysis layer (§2). No new physics. | Bitwise parity + all lockstep suites green | 2–3 days |
| **S2** | `ShellSection`, membrane-only quad with drilling DOF. Kernel layout/perf spike. | Membrane benchmarks, coarse-mesh wall test | 3–4 days |
| **S3** | MITC4 plate → full `ShellElement`; local frame, material axes, warping warning. | Plate + shell benchmarks | 3–5 days |
| **S4** | Area loads, shell self-weight, resultants, `section_cut`, pier/spandrel helper, cases/combos. | Equilibrium tests, hand-checked cut forces | 3–4 days |
| **S5a** | Wall mesher (rectangle + rectangular openings, structured, hard points at storey levels and frame nodes, edge-node sets returned); a **2-node link/spring element** for panel joints. | Mixed wall-frame building example (README-tested) | 2–3 days |
| **S5b** | Rigid diaphragm constraints: per-element transformation K̂ₑ = Tₑᵀ Kₑ Tₑ (DERIVATION §5), constrained-DOF class in `DofPartition`, load transformation, slave recovery, reactions. | Rigid-diaphragm frame ≡ very-stiff-slab model; torsional load-sharing hand check | 3–4 days |
| **S6** | Pure path / AD: shells (and diaphragm Tₑ) in `ModelState` and `stiffness_entries`; gradient tests; AsapOptim `ThicknessVariable` + shell-aware `OptParams`. | AD checks; AsapOptim suite | 2–3 days |
| — | **Release Asap 1.3.0** (S1–S6) | registry | |
| **S7** | Generic shells: 3-node element, mesh import (Gmsh/Meshes.jl via extension, no hard deps), warped-quad correction, MITC4+ if needed. | Shell benchmark set on triangles and mixed meshes | later |

## 7. Companion features

- **Rigid diaphragm constraints (IN 1.3, phase S5b).** Multipoint constraints by
  transformation. Shell floors already give a *semi*-rigid diaphragm; the rigid constraint is
  the cheap, conventional option. Asap has no constraint mechanism today.
- **Precast joint modelling**: a continuous shell mesh across precast panel joints overestimates
  wall stiffness. Real precast shear walls are governed by their horizontal/vertical connections.
  The S5 link element is the minimum; a joint "interface" helper (duplicate nodes along a line,
  connected by springs) is the convenient version.
- Mass matrices (dynamics hook) and geometric stiffness (P-Δ) for shells ride with the deferred
  nonlinear/dynamics phase.

## 8. Decisions

Resolved 2026-11-23; see the Decisions section at the top.

**Downstream note (S6):** AsapOptim builds its own `i1`/`i2`/`Cinc` from
`el.nodeStart`/`el.nodeEnd` over *all* `model.elements` (`parameters.jl`). A model containing
shells must either be supported there (shell groups excluded from truss/frame batched paths) or
rejected with a clear error. That decision belongs to S6, and AsapOptim's suite joins the gate.

## References (to be pulled in S0)

- Ibrahimbegović A., Taylor R.L., Wilson E.L. (1990). A robust quadrilateral membrane finite
  element with drilling degrees of freedom. *IJNME* 30(3):445–457. doi:10.1002/nme.1620300305
- Hughes T.J.R., Brezzi F. (1989). On drilling degrees of freedom. *CMAME* 72:105–121.
- Dvorkin E.N., Bathe K.J. (1984). A continuum mechanics based four-node shell element for general
  nonlinear analysis. *Engineering Computations* 1:77–88.
- Bathe K.J., Dvorkin E.N. (1985). A four-node plate bending element based on Mindlin/Reissner
  plate theory and a mixed interpolation. *IJNME* 21:367–383.
- Batoz J.L., Tahar M.B. (1982). Evaluation of a new quadrilateral thin plate bending element.
  *IJNME* 18:1655–1677.
- Batoz J.L., Bathe K.J., Ho L.W. (1980). A study of three-node triangular plate bending elements.
  *IJNME* 15:1771–1812.
- Allman D.J. (1984). A compatible triangular element including vertex rotations for plane
  elasticity analysis. *Computers & Structures* 19:1–8.
- MacNeal R.H., Harder R.L. (1985). A proposed standard set of problems to test finite element
  accuracy. *Finite Elements in Analysis and Design* 1:3–20.
- Ko Y., Lee P.S., Bathe K.J. (2017). A new MITC4+ shell element. *Computers & Structures*
  182:404–418.
- OpenSees `ShellDKGQ` (GQ12 membrane with drilling + DKQ plate; Wang, Lu, Xie, Cen, Gu) and
  Lu X.Z. et al. (2015), *FEAD* 98:14–25 (nonlinear shear-wall shell in OpenSees).
- Timoshenko S., Woinowsky-Krieger S. *Theory of Plates and Shells* (plate reference values).

(Citations from memory except Ibrahimbegović 1990 and the OpenSees DKGQ docs, which were checked
2026-11-23. All get verified against the PDFs in S0.)

### S6d — reverse-mode performance (overnight, 2026-11-23/24)

Keith: "carry forward with any cleanup of the code, any performance gains". Profiled first
(`examples/ad/reverse_scaling.jl`):
- Zygote frame gradients cost ~150× the primal: per element, the stiffness pullback was 150 µs
  (kernel 0.8 µs) and the element-load pullback 256 µs (kernel 0.1 µs) — Zygote tracing 12×12 SMatrix
  kernels. Shell gradient per element grew with model size (650 → 1266 µs, 96 → 1536 el): O(n²) from
  captured-array indexing in the per-element closures.
- KEPT `af97018` — `AsapFrameADExt`: local ForwardDiff rules on plain-data frame kernels (bitwise
  primal). Frames 62.8 → 16.7 ms (120 el), 436 → 256 ms (665 el).
- REVERTED `32f23aa` → `f8aa055` — gather + `map` over gathered arrays fixed the O(n²) for Zygote
  (shells 2.5×, frames to 11–13× with the rules) but Enzyme mishandles `map` over `Vector{Any}`:
  reverse throws, FORWARD SILENTLY WRONG. Caught only by the Enzyme matrix rerun.
- Verified on the final code: Asap suite, lockstep, Enzyme matrix 1.13 + 1.12, regression vs v1.2.0
  (primal + ForwardDiff gradients bitwise, Zygote ≤ 6e-15).
- DONE (2026-11-24, Keith: "get the safer shell implementation for Enzyme and Zygote performance
  done before release") — Enzyme-safe linear-cost Zygote: RuleConfig rrules in AsapChainRulesExt on
  the per-element maps (`_group_entries_all` for frame groups, `_group_entries_nnode_all` for shell
  groups, `_frame_fef`, `_shell_fef`, `fef_local`). Primal code unchanged (one helper,
  `_frame_entries`, factored out); the pullback calls `rrule_via_ad` per element on its own gathered
  arguments and accumulates ONE dX / sections / ends cotangent. Enzyme and Mooncake import only
  config-free kernel rules, so they never see these. `_natural` converts structural
  `Tangent{SMatrix}` fields back to arrays.
- DONE — shell kernel pullbacks: separate ForwardDiff passes per input group (membrane rigidities,
  plate rigidities, load, positions). Only the position group carries dual geometry: stiffness
  pullback 200 → 125 µs, compile ~85 → 43 s. `_shell_rigidities(sec)` (one rule for `ShellSection`:
  one width-4 dual pass over E, G, ν, t, instead of a 35 µs Zygote trace of four accessors).
  `_shell_fef_plain(A, qg, ql, X…)` takes the traction as global + local parts (the kernel forms
  qg + Rᵀql, bitwise the old traction), so Zygote no longer traces `shell_frame` per loaded element.
- FIXED (pre-existing) — Zygote through `ShellRigiditySection(sec; m11 …)`: the same-type
  `SMatrix{3,3,T}(::SMatrix)` conversion has no clean pullback; skipped via `_smat`.
