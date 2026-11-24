# Shell walkthroughs

Shift+enter friendly scripts: each `begin … end` block is one step, prose
in `#= … =#` blocks explains it, and every file also runs top to bottom.
Units: kN, m.

| File | Shows |
|---|---|
| `01_wall_basics.jl` | A shear wall with a door and a window, **standalone**: meshing with openings, self-weight + storey loads, drift, membrane-force contours, compression struts (principal trajectories), pier forces by `section_cut`, mesh convergence |
| `02_slab_plate.jl` | A slab, **standalone** (plate bending): simply supported slab vs Timoshenko (thin vs real thickness), a flat plate on point supports, deflection and moment contours |
| `03_frame_wall_building.jl` | Shells **with frames**: perimeter frames + a coupled shear wall (frame coupling beams in the wall's plane, through the drilling DOF), rigid diaphragms, load cases and combinations, drift, wall-vs-frame base-shear split, pier and coupling-beam forces, envelopes, 3D figures |
| `04_precast_joints.jl` | Stacked precast panels joined by `LinkElement`s: drift vs joint shear stiffness, joint axial-force distribution (tie demand) |
| `05_thickness_design.jl` | Design sensitivities (Zygote, implicit `solution_tangents`, finite differences) and an optimality-criteria thickness redistribution with AsapOptim's `ThicknessVariable` |

`../ad/ad_backends_shells.jl` compares Zygote, Enzyme (reverse) and Mooncake through shells on the same
models (all agree to 1e-16; Enzyme fastest).

`figures.jl` holds the figure generators (`fig_wall`, `fig_principal`,
`fig_slab`, `fig_3d`, `fig_bars`, `fig_series`) — plain functions of Asap's
`ShellGeo`/`WallMesh` — and every walkthrough saves its figures to
`figures/`.

## Running

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples examples/shells/01_wall_basics.jl
```

The shared `examples/` environment uses Asap and AsapOptim from this checkout (`[sources]`
paths), CairoMakie, and the kjlMakie theme. The first run compiles the
plotting stack (a few minutes); later runs are fast.
