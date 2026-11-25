# Load cases, combinations, envelopes

```@setup cases
using Asap, LinearAlgebra
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

## On the differentiable path

The pure solve ([`solve`](@ref)`(model, state)`, the AD entry point) sees one load vector. For several
cases at one design, build a [`CaseLoads`](@ref) table once with [`case_loads`](@ref), then solve every
case at any design state with one factorization:

```@example cases
cl = case_loads(model)                       # every case; or case_loads(model; cases = [:dead, :live])
state = extract_state(model)
U = solve(model, state, cl)                  # n_dof × 3, one column per case
compliance(model, state, cl)                 # one value per case
```

Each case's element loads are re-evaluated at the design (`loads = :design`, the default): a line load
follows its member's length, self-weight the section, a shell pressure the element's area and
orientation. `loads = :frozen` keeps them at the model's own state.

A combination is a linear combination of the columns, with weights from
[`combination_weights`](@ref). Its compliance is the work of the combined load on the combined
displacement, not the weighted sum of the case compliances (compliance is quadratic in the loads):

```@example cases
λ = combination_weights(cl, strength)        # [1.2, 1.6, 0.5]
u = U * λ                                    # ≡ combine(cr, strength).u
free = model.cache.partition.free
F = (cl.P - case_Pf(cl, state))[free, :]     # per-case loads at the state
c = dot(F * λ, u[free])                      # ≡ combine(cr, strength).compliance
```

All of it differentiates with Zygote, ForwardDiff, Enzyme and Mooncake, with respect to positions, sections, shell thicknesses and end springs. The per-case frame fixed-end
forces for force recovery are [`fef_local`](@ref)`(cl, state)`: column `j` belongs to case
`cl.names[cl.fcase[j]]` and element `cl.fef.elem[j]`.
