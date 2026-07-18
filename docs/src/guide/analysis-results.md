# Analysis and results

```@setup analysis
using Asap
```

Build a model, solve it, query the results:

```@example analysis
steel  = Material(200e6, 80.0, 0.3)
column = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)

n1 = Node([0.0, 0.0, 0.0], :fixed)
n2 = Node([0.0, 0.0, 3.0], :free)
el = FrameElement(n1, n2, column)

model = Model([n1, n2], [el], [NodeForce(n2, [10.0, 0.0, 0.0])])
solve!(model)                     # process! runs automatically the first time
```

Results are read through accessor functions — they never live on nodes or elements:

```@example analysis
res = model.results
displacement(res, n2)             # SVector{6} at a node
```

```@example analysis
reaction(res, n1)                 # support reactions (incl. moment reactions)
```

```@example analysis
element_forces(res, el)           # local end-force 12-vector
```

```@example analysis
axial_force(res, el)              # tension-positive scalar
```

```@example analysis
res.compliance                    # external work Fᵀu
```

## Repeated solves

Repeated solves — geometry or section *values* changed, topology unchanged — reuse the frozen sparsity pattern and refactorize numerically: just call `solve!(model)` again. After topology changes (elements added, releases toggled, supports changed): `solve!(model; reprocess = true)`.
