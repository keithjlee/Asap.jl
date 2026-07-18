# Nodes and supports

```@setup nodes
using Asap
```

Every node has six DOFs — three translations, three rotations — with a fixity per DOF (`true` = free). Common supports have symbols; anything else is an explicit 6-vector:

```@example nodes
base   = Node([0.0, 0.0, 0.0], :fixed)      # clamped
pin    = Node([5.0, 0.0, 0.0], :pinned)     # translations fixed, rotations free
roller = Node([10.0, 0.0, 0.0], [true, true, false, true, true, true])  # z-roller
fixnode!(roller, :zfixed)                   # change a support later
roller
```

The full set of support symbols lives in [`FIXITIES`](@ref).

`planarize!(model)` constrains a model to 2D behavior (fixes out-of-plane DOFs and zeroes roll angles for the `:XY` plane).

## Inactive DOFs

Rotational DOFs that nothing stiffens — a node connected only to truss members, or only to fully released member ends — are classified **inactive** and excluded from the solve automatically. No singular stiffness matrices, no manual bookkeeping: an all-truss model solves a translations-only system by construction.
