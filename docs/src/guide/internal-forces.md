# Internal forces and deflected shapes

```@setup recovery
using Asap
steel  = Material(200e6, 77e6, 8.0, 0.3)
wshape = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)
rb1 = Node([0.0, 0.0, 0.0], :fixed)
rb2 = Node([6.0, 0.0, 0.0], :free)
beam = FrameElement(rb1, rb2, wshape, :beam)
model = Model([rb1, rb2], AbstractElement{Float64}[beam],
    AbstractLoad{Float64}[LineLoad(beam, [0.0, 0.0, -2.0]),
        PointLoad(beam, 0.4, [0.0, 0.0, -10.0])])
solve!(model)
```

Recovery is **equilibrium-based**: starting from the exact member end actions, the applied loads are integrated analytically along the member — the recovered fields are exact for any end condition (releases, semi-rigid springs) with no per-case formulas.

The examples below use a solved cantilever `model` with a `beam` carrying a uniform [`LineLoad`](@ref) and a mid-span [`PointLoad`](@ref).

## Scalar queries

At any fraction of the length, with zero allocation:

```@example recovery
moment_z(model, beam, 0.5)     # bending about local z — pairs with Vy, sagging+ for +y loads
```

```@example recovery
shear_y(model, beam, 0.25), axial_force(model, beam, 0.0), torsion(model, beam, 0.5)
```

Names are axis-correct: `Mz` is the moment about local z and pairs with `Vy` (`dMz/dx = Vy`); `My` pairs with `Vz` (`dMy/dx = −Vz`).

## Dense sampling for plotting

Stations include every load breakpoint and **both** sides of each point action, so shear jumps render as true jumps:

```@example recovery
f = InternalForces(model, beam; resolution = 40)
f.x'
```

```@example recovery
[f.N f.Vy f.Mz f.Vz f.My f.Mx]'
```

## Deflected shapes

Exact local displacements at any fraction:

```@example recovery
st = internal_forces(model, beam)
local_displacements(st, 0.5)   # (u, v, w) in local axes
```
