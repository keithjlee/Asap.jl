# Elements

```@setup elements
using Asap
steel  = Material(200e6, 77e6, 8.0, 0.3)
wshape = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)
bar    = Section(steel, 5e-3)
```

The examples below reuse `steel`, `wshape` (a frame section), and `bar` (an axial-only section) from [Materials and sections](materials-sections.md).

## Frame elements

3D Euler-Bernoulli beam-columns carrying axial force, biaxial bending, shear, and torsion. `rollangle` is the section roll angle about the member axis (default `π/2`).

```@example elements
na = Node([0.0, 0.0, 0.0], :fixed)
nb = Node([4.0, 0.0, 0.0], :free)

girder = FrameElement(na, nb, wshape, :girder)
rolled = FrameElement(na, nb, wshape; rollangle = 0.0)
girder
```

### End releases and semi-rigid connections

Connections are **end springs** in the element's local axes. The classical releases are exact limits (`Inf` = rigid, `0` = released), and any finite value is a semi-rigid connection:

```@example elements
hinged = FrameElement(na, nb, wshape; release = :fixedfree)   # hinge at far end
# available: :fixedfixed (default), :fixedfree, :freefixed, :freefree, :joist

# semi-rigid: e.g. a bolted end plate with finite rotational stiffness
semi = FrameElement(na, nb, wshape,
    EndConditions(EndSprings(Inf, Inf, 5e4, 5e4), rigid_end()), :connection)
```

Connection stiffness is ordinary data — it can be a design variable in gradient-based optimization, and everything downstream (fixed-end forces, force recovery) handles it with no special cases.

## Truss elements — and mixing them with frames

Axial-only two-force members. They coexist freely with frame elements in one model; a truss element simply never touches rotational DOFs:

```@example elements
n1 = Node([0.0, 0.0, 0.0], :fixed)
n2 = Node([3.0, 0.0, 0.0], :free)
n3 = Node([3.0, 0.0, 3.0], :fixed)                # an anchor above

beam = FrameElement(n1, n2, wshape, :beam)
tie  = TrussElement(n2, n3, bar, :tie)

model = Model([n1, n2, n3],
    AbstractElement{Float64}[beam, tie],          # mixed vector: type it explicitly
    AbstractLoad{Float64}[NodeForce(n2, [0.0, 0.0, -100.0])])
solve!(model)
axial_force(model.results, tie)                   # tension-positive
```

## Variable elements (varying cross-section)

One user-facing member whose section varies along its length as a chain of prismatic segments — for haunched girders, stepped columns, or optimization results. Interior joints become internal DOFs: the model is never mutated, no phantom nodes appear, and you query the member as a single piece:

```@example elements
deep    = Section(steel, 2e-2, 4e-4, 1e-4, 2e-6)
mid     = Section(steel, 1.5e-2, 2e-4, 6e-5, 1e-6)
shallow = wshape

v1 = Node([0.0, 0.0, 0.0], :fixed)
v2 = Node([8.0, 0.0, 0.0], :free)
haunched = VariableElement(v1, v2,
    AbstractSection{Float64}[deep, mid, shallow],   # one section per segment
    [0.25, 0.6])                                    # interior break fractions

vmodel = Model([v1, v2], AbstractElement{Float64}[haunched],
    AbstractLoad{Float64}[LineLoad(haunched, [0.0, 0.0, -2.0])])
solve!(vmodel)
moment_z(vmodel, haunched, 0.5)                 # resolves to the right segment
```

```@example elements
element_forces(vmodel.results, haunched, 2)     # or reach a specific segment
```

`SelfWeight` on a `VariableElement` automatically varies with each segment's `ρA`:

```@example elements
push!(vmodel.loads, SelfWeight(haunched; g = [0.0, 0.0, -9.81]))
solve!(vmodel)
moment_z(vmodel, haunched, 0.5)
```
