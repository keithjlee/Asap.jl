# Force density method (form-finding)

```@setup fdm
using Asap
```

A self-contained FDM subsystem for cable/tension networks ships in the same package, built on the same three-layer architecture as the frame core: [`Network`](@ref) holds pure definition data, `process!` freezes the force-density pattern once per topology, and repeated solves are a pure scatter of the current `q` values plus a numeric-only refactorization — so form-finding sweeps are cheap by construction. `solve!` updates the node positions in place (the geometry *is* the result); forces and reactions are queried from `network.results`:

```@example fdm
# a 2-cable chain: anchors at the ends, load at the middle
a = FDMnode([0.0, 0.0, 0.0], false)          # false = fully fixed (anchor)
m = FDMnode([1.0, 0.0, 0.0], true, :mid)     # true  = fully free
b = FDMnode([2.0, 0.0, 0.0], false)

cables = [FDMelement(a, m, 2.0), FDMelement(m, b, 2.0)]   # q = 2.0 [force/length]
net = Network([a, m, b], cables, [FDMload(m, [0.0, 0.0, -1.0])])

solve!(net)
m.position                            # the found geometry
```

```@example fdm
member_force(net.results, cables[1])  # cable tension q·L
```

```@example fdm
reaction(net.results, a)              # anchor force, per fixed axis
```

Form-finding sweep: mutate `q` and re-solve (numeric refactorization only):

```@example fdm
update_q!(net, 4.0)                   # doubles q — halves the sag
m.position
```

## Per-axis fixity

Fixity is **per-axis** (`true` = free, like `Node`): prescribe the plan and let the height find equilibrium — natural for vault form-finding. FDM equilibrium is separable per coordinate, so each axis solves against its own free/fixed partition:

```@example fdm
crown = FDMnode([1.0, 1.0, 0.0], Bool[false, false, true])   # plan fixed, z free
```

The `solver` keyword works here too (any LinearSolve algorithm once loaded). [`to_network`](@ref)`(model)` converts a solved truss model into an equivalent FDM network, and [`to_truss`](@ref)`(network, section)` converts back. Differentiable force-density optimization (`QVariable`, `solve_network`) lives in [AsapOptim](https://github.com/keithjlee/AsapOptim).
