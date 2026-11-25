# Unreleased (after v1.4.0)

## Added
- Load cases on the differentiable path. `case_loads(model; cases)` builds a `CaseLoads` table once;
  `solve(model, state, cl)` returns every case's displacements (`n_dof × n_cases`) from ONE stiffness
  factorization, `compliance(model, state, cl)` the per-case compliances. Every case's element loads
  (frame line/point loads and self-weight, quad and triangle shell surface loads) are re-evaluated at
  the design (`loads = :design`) or kept at the model's state (`:frozen`); rigid diaphragms are
  supported. `combination_weights(cl, combo)` turns a `LoadCombination` into the column weights λ
  (`U * λ`); `case_Pf(cl, state)` gives the per-case fixed-end forces and `fef_local(cl, state)` the
  per-load local ones for force recovery. Differentiable with Zygote (linear-cost rules) and
  ForwardDiff.
- `mesh_polygon(…; segments)`: interior lines the element edges must follow, given as point pairs
  `(a, b)` (section cuts through the material, interior load or bearing lines). Their nodes come from
  `nodes_on_segment`, so `section_cut` integrates the exact force across any chosen line. A segment
  may end on the outline, a hole or another segment; a segment leaving the material, or crossing an
  edge or another segment, is refused.
- Walkthrough `examples/shells/07_dapped_end.jl`: a dapped-end beam. Bearing pads as pressures with
  rigid-body supports only, exact interface cut forces, the nib tie and the hanger as integrated cut
  forces, the crack direction at the re-entrant corner, and a mesh study.
