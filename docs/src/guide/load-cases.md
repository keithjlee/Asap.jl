# Load cases, combinations, envelopes

```@setup cases
using Asap
steel  = Material(200e6, 77e6, 8.0, 0.3)
wshape = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)
```

Tag loads with cases; solve every case against **one** factorization; get any factored combination by superposition — no re-solves, ever.

The example is a small portal frame with dead, live, and wind cases:

```@example cases
n1 = Node([0.0, 0.0, 0.0], :fixed)
n2 = Node([0.0, 0.0, 3.0], :free)
n3 = Node([5.0, 0.0, 3.0], :free)
n4 = Node([5.0, 0.0, 0.0], :fixed)
beam = FrameElement(n2, n3, wshape, :beam)
elements = AbstractElement{Float64}[
    FrameElement(n1, n2, wshape), beam, FrameElement(n4, n3, wshape)]

loads = AbstractLoad{Float64}[
    LineLoad(beam, [0.0, 0.0, -2.0]; case = :dead),
    PointLoad(beam, 0.4, [0.0, 0.0, -10.0]; case = :live),
    NodeForce(n2, [5.0, 0.0, 0.0]; case = :wind),
]
model = Model([n1, n2, n3, n4], elements, loads)

cr = solve_cases!(model)                    # one assembly, one factorization
```

Combinations are exact, by superposition:

```@example cases
strength = LoadCombination(:LRFD, [:dead => 1.2, :live => 1.6, :wind => 0.5])
res = combine(cr, strength)
displacement(res, n2)
```

Station-wise min/max over a combination set — what design checks consume:

```@example cases
env = envelope(model, beam, cr,
    [strength, LoadCombination(:service, [:dead => 1.0, :live => 1.0])])
env.x, env.lo, env.hi                       # rows: N, Vy, Mz, Vz, My, Mx
```
