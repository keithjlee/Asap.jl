# Loads and springs

```@setup loads
using Asap
steel  = Material(200e6, 77e6, 8.0, 0.3)
wshape = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)
b1 = Node([0.0, 0.0, 0.0], :fixed)
n2 = Node([6.0, 0.0, 0.0], :free)
beam = FrameElement(b1, n2, wshape, :beam)
```

## Loads

All loads take an `id` and a `case` tag (default `:LC1`). Given a node `n2` and a frame element `beam`:

```@example loads
NodeForce(n2, [0.0, 0.0, -50.0]; case = :live)
NodeMoment(n2, [0.0, 1e3, 0.0])

LineLoad(beam, [0.0, 0.0, -2.0])                            # uniform, full span
TrapezoidLoad(beam, 0.2, 0.8, 1.0, 3.0, [0.0, 0.0, -1.0])   # partial-span, varying
DistributedLoad(beam, [0.0, 0.3, 0.5], [0.0, 4.0, 0.0],     # arbitrary piecewise-linear
    [0.0, 0.0, -1.0])
PointLoad(beam, 0.4, [0.0, 0.0, -10.0])                     # at 40% of the span
PointMoment(beam, 0.5, [0.0, 0.0, 800.0])                   # concentrated moment

SelfWeight(beam; g = [0.0, 0.0, -9.81])                     # from ρA(section)
```

Every distributed shape lowers to one canonical piecewise-linear type with one exact integration (3-point Gauss against the element shape functions — exact, not approximate, for these loads). Adding a custom element load type means implementing a single method, [`fixed_end_forces`](@ref).

## Spring supports

Elastic supports are **applicative**: a [`NodalSpring`](@ref) references its node (the way a load does) and lives on the model — nodes don't know about their springs, and several springs on one node add up.

```@example loads
sb  = Node([0.0, 0.0, 0.0], :fixed)
tip = Node([0.0, 0.0, 3.0], :free)

soil = NodalSpring(tip, [0.0, 0.0, 5e4, 0.0, 0.0, 0.0], :soil)   # vertical only
pad  = NodalSpring(tip, 1e5)                                     # uniform translational

smodel = Model([sb, tip], AbstractElement{Float64}[FrameElement(sb, tip, wshape)],
    AbstractLoad{Float64}[NodeForce(tip, [1.0, 0.0, -10.0])]; springs = [soil])
solve!(smodel)
displacement(smodel.results, tip)
```

Spring reactions are recovered as `−k·u` in post-processing.
