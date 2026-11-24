# Shells, walls and diaphragms

Shell sections and elements, their kernels, surface loads, resultants and section cuts, the wall mesher, link elements and rigid diaphragms.

```@autodocs
Modules = [Asap]
Pages = ["Materials/shell_sections.jl", "Elements/shell.jl", "Elements/link.jl",
         "Elements/kernels/shell.jl", "Elements/kernels/shell_loads_recovery.jl",
         "Loads/shell_loads.jl", "Results/shell_results.jl", "Generation/WallMesh.jl",
         "Model/diaphragm.jl", "Analysis/constraints.jl", "Geometry/ShellGeometry.jl"]
Private = false
```

## Internals

```@docs
Asap.ConstraintMap
Asap.mitc4_plate_stiffness
Asap.shell_load_vector
Asap.shell_consistent_load
Asap.shell_resultants_local
Asap.build_constraint
```
