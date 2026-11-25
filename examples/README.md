# Asap walkthroughs

Shift+enter friendly scripts. Each `begin … end` block is one step, the prose in `#= … =#`
blocks explains it, and every file also runs top to bottom. All files share one environment
(this folder's `Project.toml`) and one set of figure generators (`figures.jl`). Figures are
saved to `figures/` inside each subfolder.

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples examples/analysis/01_truss_2d.jl
```

Units in every walkthrough: kN, m. Asap itself is units-agnostic.

## Analysis — `analysis/`

| File | Scenario |
|---|---|
| `01_truss_2d.jl` | Planar Pratt truss, built by hand: supports, `planarize!`, axial forces (statics-checked), reactions, deflected shape; the `Warren2D` generator; a depth study |
| `02_spaceframe_3d.jl` | 3D double-layer space truss (`SpaceFrame`): self-weight lumped to nodes (trusses take no element loads), reactions = applied load, force distribution, support layouts |
| `03_frame_2d.jl` | 2D portal frame: line, point and lateral loads, section orientation (roll angle), member N/V/M diagrams, exact deflected shape, end releases (pinned girder → wL²/8) |
| `04_frame_3d_building.jl` | 3D building frame (`Frame`): load cases D/L/W solved against one factorization, storey drift, combinations, a member moment envelope |
| `05_joints_springs_variable.jl` | Semi-rigid connections (`EndSprings` sweep between the pinned and rigid limits), spring supports (`NodalSpring` footing; matches the closed form), a haunched girder (`VariableElement`) |

## Shells, walls and diaphragms — `shells/` (see `shells/README.md`)

1. A shear wall with openings.
2. Slabs and plates.
3. A frame plus coupled shear wall building with rigid diaphragms.
4. Precast panel joints.
5. Wall-thickness sensitivities and optimization.
6. A disturbed region: a deep transfer panel with openings, triangulated (`mesh_polygon`).

## Optimization — `optimization/` (with AsapOptim)

| File | Scenario |
|---|---|
| `01_truss_size.jl` | SIZE: member areas (`AreaVariable`), minimum compliance at fixed volume by optimality criteria; ends fully stressed (uniform stress in every active member) |
| `02_truss_shape.jl` | SHAPE vs SIZE: symmetric top-chord heights (`SpatialVariable` + `CoupledVariable`) against areas, compared on C·V (independent of the amount of steel) |
| `03_frame_joints.jl` | Connection stiffness (`JointVariable`): minimum total joint stiffness at drift = H/400, gradients by implicit differentiation (`solution_tangents`) |

## Automatic differentiation — `ad/`

- `enzyme_matrix.jl` + `run_enzyme_matrix.sh`: Enzyme (reverse and forward) against Zygote
  (or finite differences) on every scenario class — trusses, frames with element loads,
  semi-rigid joints, springs, variable elements, shells, rigid diaphragms, a large space truss.
  Results are in `ad/ENZYME_STATUS.md`.
- `ad_backends_shells.jl`: Zygote vs Enzyme vs Mooncake timings through shell models.
- `reverse_scaling.jl`: reverse-mode (Zygote) cost against model size for frames and shells, with
  the v1.2 → 1.3 numbers for frames.

## Reports — `reports/`

One PDF report per walkthrough (`A1_truss_2d.pdf` … `D2_ad_backends_shells.pdf`) and `main.pdf`, all
of them in one document. Each walkthrough ends with a `report_data("<id>"; …)` block (`figures.jl`)
that records the numbers its report prints, so no number in a report is typed by hand.

