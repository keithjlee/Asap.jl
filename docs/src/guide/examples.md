# Walkthroughs

The repository's `examples/` folder has shift+enter walkthroughs for every class of structure
Asap handles. Each is a script of `begin … end` steps with prose between them, and each also
runs top to bottom. They share one environment and a set of figure generators.

- **Analysis:**
  - planar and 3D trusses;
  - 2D portal frames (loads, releases, member diagrams, exact deflected shapes);
  - 3D building frames (load cases, combinations, envelopes);
  - semi-rigid joints, spring supports, variable-section members.
- **Shells, walls and diaphragms:**
  - shear walls with openings, slabs and plates;
  - a frame plus coupled-wall building with rigid diaphragms;
  - precast joints;
  - wall-thickness design.
- **Optimization (AsapOptim):** truss size, shape vs size, connection stiffness.
- **Automatic differentiation:** an Enzyme verification matrix across all scenario classes,
  and backend timings.

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples examples/analysis/01_truss_2d.jl
```

See `examples/README.md` for the full index.
Each walkthrough also has a Typst report in `examples/reports/`, which prints the numbers the
walkthrough writes to `examples/reports/data/` when it runs.
