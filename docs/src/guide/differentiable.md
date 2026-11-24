# Differentiable analysis

```@setup diff
using Asap
steel  = Material(200e6, 77e6, 8.0, 0.3)
wshape = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)
n1 = Node([0.0, 0.0, 0.0], :fixed)
n2 = Node([0.0, 0.0, 3.0], :free)
n3 = Node([5.0, 0.0, 3.0], :free)
n4 = Node([5.0, 0.0, 0.0], :fixed)
beam = FrameElement(n2, n3, wshape, :beam)
elements = AbstractElement{Float64}[
    FrameElement(n1, n2, wshape), beam, FrameElement(n4, n3, wshape)]
loads = AbstractLoad{Float64}[
    LineLoad(beam, [0.0, 0.0, -2.0]),
    NodeForce(n2, [5.0, 0.0, 0.0]),
]
model = Model([n1, n2, n3, n4], elements, loads)
```

The entire pipeline has a **pure functional path**: extract a [`ModelState`](@ref) (positions and sections as plain data), evaluate, differentiate. Loading Zygote (or anything ChainRules-aware) activates Asap's rule extension automatically; a handful of small rules cover the sparse construction and the linear solve, and everything else differentiates natively.

The example differentiates the compliance of a small portal frame `model` with respect to every node coordinate:

```@example diff
using Zygote

solve!(model)                               # builds the analysis cache
state = extract_state(model)

# gradient of compliance w.r.t. every node coordinate — a 3 × n sensitivity field:
g = Zygote.gradient(state.X) do X
    compliance(model, ModelState{Float64}(X, state.sections))
end[1]
```

Gradients flow with respect to node positions, section properties, and semi-rigid connection stiffnesses.

For fast Zygote gradients, load ForwardDiff as well. The frame and shell kernels then get reverse
rules that evaluate each element's derivative with ForwardDiff (`AsapFrameADExt`,
`AsapShellADExt`) instead of tracing the small static-matrix kernels scalar by scalar (frames: about
4× faster at 120 elements). Without ForwardDiff, the gradients are the same, just slower. See
`examples/ad/reverse_scaling.jl` for timings against model size.

## Element loads follow the design

Element-load fixed-end forces are re-evaluated at the state ([`assemble_Pf`](@ref)), so the load path is part of every gradient. In the portal above, moving a knee changes the beam's length and with it the total line load, and the gradient includes that term. The defaults under design change:

| Load | Under a design change |
|---|---|
| `LineLoad` / `DistributedLoad` | per-unit-length intensity over fractions of the member: the total is recomputed from the current length |
| `PointLoad` / `PointMoment` | stays at its fraction of the member, riding with it |
| `SelfWeight` | follows the current ρA and length |
| `:local` directions | rotate with the member (`:global` ones do not) |
| semi-rigid / released ends | condensed end moments follow the end springs and EI |

Pass `loads = :frozen` to `solve`/`compliance` to keep the fixed-end forces at the model's stored design (the behavior before v1.2), for studies where the load is fixed in space. Nodal loads (`NodeForce`, `NodeMoment`) are constant either way. For design-variable bookkeeping (areas, coupled symmetric geometry, bounds) and optimization-ready objectives, use [AsapOptim](https://github.com/keithjlee/AsapOptim) — a thin layer over this path, verified against the original published implementation to 13 digits (see its `docs/AD_VERIFICATION_AND_BENCHMARKS.md`).

Forward-mode AD (ForwardDiff, Enzyme forward) and Mooncake are supported through the same path via their own package extensions — load the AD package and differentiate; no further setup.
