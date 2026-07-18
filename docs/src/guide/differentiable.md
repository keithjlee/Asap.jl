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

Gradients flow with respect to node positions, section properties, and semi-rigid connection stiffnesses. For design-variable bookkeeping (areas, coupled symmetric geometry, bounds) and optimization-ready objectives, use [AsapOptim](https://github.com/keithjlee/AsapOptim) — a thin layer over this path, verified against the original published implementation to 13 digits (see its `docs/AD_VERIFICATION_AND_BENCHMARKS.md`).

Forward-mode AD (ForwardDiff, Enzyme forward) and Mooncake are supported through the same path via their own package extensions — load the AD package and differentiate; no further setup.
