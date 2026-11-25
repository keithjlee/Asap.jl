# S7: triangles and general meshes — plan and tracker

> **Status: APPROVED 2026-11-24 (Keith). In progress.** Branch `triangles` (off `main` @ `8a49467`,
> Asap 1.3.0). Target: Asap **1.4.0** (additive). Follows `PLAN.md` (S0–S6, released as 1.3).
> Formulations and benchmark values: `DERIVATION.md` §7 onward (added as each element lands).

## Why now

Strut-and-tie modelling of disturbed regions (panels, deep beams, dapped ends, corbels, walls with
openings) needs a **plane-stress solve of an arbitrary polygon with polygon or circular holes**:
principal stress magnitudes and directions (load-path tracing, stress-field underlays), support
reactions per bearing, and envelopes over load cases. `mesh_wall` covers rectangles with
rectangular openings only. The general case needs triangles and a general mesher. Higher-order
elements, adaptive refinement and nonlinear analysis are not required.

## Decisions (Keith, 2026-11-24)

1. **Membrane: ANDES-OPT** (Felippa 2003, CMAME 192:2125). The optimal 3-node drilling triangle:
   exact in-plane pure bending of rectangular mesh units, closed form, no spurious modes. Same
   element as OpenSees `ShellANDeS` (independent cross-check).
   - `K = V⁻¹LELᵀ + (3/4)β₀ T̃ᵀ_θu K_θ T̃_θu` (Felippa eqs 21–32), α_b = 3/2, β₁…β₉ =
     (1, 2, 1, 0, 1, −1, −1, −1, −2) (eq 42).
   - Isotropic β₀ = max(½(1−4ν²), 0.01). General (cracked, orthotropic `ShellRigiditySection`):
     β₀ = 2/Ῡ − 3/2 with Ῡ the 2π-average of E₁₁φ·C₁₁φ (eq 44), clamped at 0.01. The average of a
     degree-8 trigonometric polynomial is exact with 16 equally spaced samples: closed form, AD-safe.
   - With h = 1 and E = A (membrane rigidity, already thickness-integrated).
2. **Plate: MITC3+** (Lee, Lee & Bathe 2014, Comput. Struct. 138:12). Reissner–Mindlin like MITC4,
   so mixed quad/triangle meshes use one plate theory and triangles report Q.
   - Rotations enriched by a cubic bubble f₄ = 27rs(1−r−s) (2 internal DOFs, condensed);
     f_i = h_i − f₄/3.
   - Assumed transverse shear (eq 17) tied at A (1/6, 2/3), B (2/3, 1/6), C (1/6, 1/6) and
     D (1/3+d, 1/3−2d), E (1/3−2d, 1/3+d), F (1/3+d, 1/3+d), d = 1/10 000.
   - Check: single-element eigenvalues of the uncondensed plate block against Table 3.
3. **Mesher: DelaunayTriangulation.jl as a weak-dependency extension** (pure Julia, constrained
   Delaunay with holes, quality refinement). Plus, in core: an array constructor (points + cells,
   triangles and quads) and a dependency-free Gmsh `.msh` (v4 ASCII) reader.

## Progress tracker

| Phase | Content | Gate | Status |
|---|---|---|---|
| S7a | N-generic shell infrastructure: `ShellElement{T,3}`, `shell_frame` for 3 nodes, FEF table node rows per element, group rules for any N, results/geo/envelope dispatch, AsapOptim shell paths | quad models bitwise identical; lockstep | ✅ in-place + pure primal (separate triangle load table `tloads…`, `_assemble_Pf_tri`); group AD rules → S7g; AsapOptim → S7g |
| S7b | ANDES-OPT membrane kernel (9×9) | patch tests; bending optimality r = 1 on rectangle pairs; Cook's membrane; rank and rigid modes | ✅ `d9c5761` (Felippa example to 5 digits; bending exact) |
| S7c | MITC3+ plate kernel (9×9 after condensation) | Table 3 eigenvalues; patch tests; clamped plate convergence; isotropy (node-order invariance) | ✅ `d9c5761` (Table 3 to 5 digits) |
| S7d | Triangle shell: 18×18 assembly, consistent loads, resultants (N, M, Q), `section_cut`, `ShellGeo`, `shell_envelope` | Scordelis–Lo and pinched cylinder on triangles; mixed quad/triangle meshes; equilibrium | ✅ Cook 1%, plates 2–3%, Scordelis–Lo 3% at 16×16; irregular patch exact; mixed mesh equilibrium + pure-path parity; `shell_faces` |
| S7e | Meshing: `shell_mesh(points, cells, section)`, `read_gmsh`, `AsapDelaunayExt` → `mesh_polygon(outline; holes, max_area, min_angle)` returning boundary edge sets | polygon with circular and L-shaped holes; mesh quality; node reuse | ✅ `read_gmsh` (MSH 4.1 ASCII, physical groups, planar or spatial), `seed` for reproducible meshes; node reuse open (`AbstractShellMesh` ⊃ `WallMesh`, `ShellMesh`) |
| S7f | D-region workflow: loads and bearings along boundary segments, patch loads, reactions per bearing, principal fields | deep beam and dapped end against hand statics; reactions sum to loads | ⏳ `nodes_on_segment`, segment `edge_loads`, `elements_in`, `total_reaction`, `fixnode!(n, ::Vector{Bool})`; deep beam statics exact. Open: dapped-end example, trajectory helper? |
| S7g | AD and optimization: plain kernels + ForwardDiff rules (width 4), Enzyme/Mooncake bridges, `ThicknessVariable` on triangles, Enzyme matrix triangle case | gradients vs FD; Enzyme matrix 1.13 + 1.12; Julia 1.10 compile time | ⏳ rules + bridges + group rules for any N; AsapOptim (branch `triangles`): mixed meshes, `shell_volume`, implicit tangents; Enzyme case `shells_tri` (rev 2e-16, fwd 3e-15); 1.10 first pullback 8.6 s (quad 141 s). Full matrix rerun pending |
| S7h | Walkthrough `examples/shells/06_d_region.jl` + PDF report; docs; release notes 1.4.0 | strict docs build; lockstep | ⏳ walkthrough written (numerics verified: exact statics, reactions converge 1298→1302→1300 kN, singular corner peaks grow); figures by the figures session (`fig_principal`/`fig_wall` for `AbstractShellMesh`); report S6 by the figures session |
| later | Warped-quad correction; MITC4+ (distorted quads) | | ☐ |

## Rules carried over from 1.3

- Every pure-path change must pass the Enzyme matrix (`examples/ad/run_enzyme_matrix.sh`, 1.13
  and 1.12) before commit. Zygote-only speedups go in `RuleConfig` rrules; the primal Enzyme sees
  stays byte-for-byte.
- ForwardDiff kernel rules use dual width 4 (wide duals compile for many minutes on Julia 1.10).
- Julia 1.10 is the floor: check compile time and the suite on 1.10 for every new kernel.
- Report sources stay private (`examples/reports-src/`); only PDFs are published.
- No client or product names in this public repo.

## Handoff notes

- **Triangle edge loads need drilling moments `∓α_b·(w·n)·l²/12`, α_b = 3/2** (the ANDES force
  lumping). Verified: only this sign and scale reproduce a uniform stress exactly (6e-14; other
  signs/scales off by 1–6 in N). Same rule as the ITW quad (α = 1).
- **A held edge of a drilling membrane must fix ψ** (as for quads): with ψ free the missing corner
  drilling moments ruin a uniform state even under exact loads.
- The quad pure-path code is untouched; triangles have their own load table fields
  (`tloads`, `telem`, `tnodes`, `tscatter`) and `_assemble_Pf_tri`. `assemble_Pf` gained one
  `isempty(t.tloads)` branch at the top, so the Enzyme matrix must be rerun.
- Zygote through triangle models works but uses no fast rules yet (the group rule is Val(4) only and
  the triangle kernels have no ForwardDiff rules): S7g.
- Test helper renamed `shell_mesh` → `grid_shell_mesh` (the export took the name).
- **Touching holes (DONE)**: `mesh_polygon` feeds DelaunayTriangulation a constrained-segment graph
  (coincident vertices merged, edges split at T-junctions) instead of boundary curves, then keeps the
  triangles whose centroid is in the material. Holes may touch the outline or each other at a vertex,
  a vertex on an edge, or along a shared edge (STM v3 openings rule: touching allowed, overlap
  refused); crossing/overlapping holes and holes outside are refused. Boundary-curve input could not
  express touching curves (DT read a curve sharing a vertex as a separate outer region).
- **Julia 1.10**: triangle kernels compile fast (primal 5.8 s, stiffness pullback 8.6 s, load
  pullback 0.2 s) — the 18×18 kernels are far lighter than the quad's 24×24.
- **Mesh reproducibility**: DelaunayTriangulation inserts points in random order; `mesh_polygon` passes a
  seeded `Xoshiro(seed)` (default 1) to `triangulate` and `refine!` (reached as
  `DelaunayTriangulation.Random`, no new Asap dependency).
- **D-region walkthrough lessons**: a pad's reaction resultant sits toward its inner edge (x = 1.075
  on a 0.4–1.2 pad), so a lever rule at pad centres is 3 % off; the local compression direction
  part-way up a pad-to-column strut is steeper than the straight resultant line (bottle-shaped
  strut, 64° vs 51°); corner peaks diverge under refinement (−3231 → −5502 kN/m).
- **Shared checkout**: another session (figures) commits to Asap `triangles`. Stage explicit paths
  only; it owns `examples/figures.jl`, walkthrough figures, report PDFs and `reports-src`.
