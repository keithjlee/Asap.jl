# Elements

Frame, truss, and variable elements; end conditions (releases and semi-rigid springs); element geometry queries and stiffness kernels.

```@autodocs
Modules = [Asap]
Pages = ["Elements/end_conditions.jl", "Elements/interface.jl", "Elements/frame.jl", "Elements/variable.jl", "Elements/kernels/transformation.jl", "Elements/kernels/stiffness.jl"]
Private = false
```

## Internal kernels

Unexported element kernels referenced by the docstrings above:

```@docs
Asap.frame_stiffness
Asap.truss_stiffness
Asap.local_stiffness
Asap.transform_to_global
Asap.segment_slots
```
