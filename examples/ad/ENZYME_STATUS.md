# Enzyme across Asap — status (2026-11-25)

**Enzyme works across Asap, in reverse and forward mode, on Julia 1.13 and 1.12.** All 20 cases
(10 scenario classes × 2 modes) pass on both Julia versions. Gradients agree with Zygote to
1e-14 or better. On `VariableElement` models, which Zygote cannot differentiate, they agree with
central finite differences to 1.6e-9, the finite-difference truncation error.

Versions: Enzyme 0.13.221, Zygote 0.7.20, DifferentiationInterface 0.7.22; Asap and AsapOptim on
branch `triangles` (Asap 1.4). Settings as everywhere in Asap: `set_runtime_activity(Reverse / Forward)` with a
`Const` function annotation, through DifferentiationInterface.

Reproduce (each case runs in its own process under `timeout 900`):

```bash
examples/ad/run_enzyme_matrix.sh                  # default julia, examples env
examples/ad/run_enzyme_matrix.sh +1.12 <project>  # another channel / environment
```

Raw logs and a JSON summary are in `results/`.

## Results

Warm medians of 3 gradient evaluations, in ms; *n* = number of design variables; error =
‖g_Enzyme − g_ref‖ / ‖g_ref‖. The objective is compliance in every case.

| Case | What is differentiated | n | Ref | Mode | Error | Enzyme 1.13 | Zygote 1.13 | Enzyme 1.12 | Zygote 1.12 |
|---|---|---|---|---|---|---|---|---|---|
| truss_area | member areas, planar truss | 14 | Zygote | rev | 8.8e-16 | 0.109 | 0.343 | 0.12 | 0.528 |
| | | | | fwd | 2.8e-15 | 0.117 | 0.334 | 0.134 | 0.524 |
| truss_shape | top-chord node heights | 3 | Zygote | rev | 6.4e-16 | 0.091 | 0.341 | 0.115 | 0.532 |
| | | | | fwd | 1.1e-15 | 0.052 | 0.313 | 0.085 | 0.667 |
| spaceframe | member areas, 8 × 8 space truss | 512 | Zygote | rev | 1.0e-14 | 0.849 | 3.9 | 0.863 | 1.3 |
| | | | | fwd | 1.5e-14 | 285 | 1.28 | 538 | 1.28 |
| frame_lineload | frame section with line load + self-weight (design-dependent FEFs) | 3 | Zygote | rev | 7.7e-15 | 0.143 | 0.569 | 0.139 | 0.743 |
| | | | | fwd | 4.1e-15 | 0.102 | 0.586 | 0.108 | 0.745 |
| frame_joints | semi-rigid connection stiffness (`JointVariable`) | 1 | Zygote | rev | 1.9e-16 | 0.48 | 0.812 | 0.487 | 1.07 |
| | | | | fwd | 3.8e-16 | 0.458 | 0.83 | 0.46 | 1.08 |
| springs | spring-support stiffness | 1 | Zygote | rev | 2.0e-16 | 0.096 | 0.302 | 0.118 | 0.393 |
| | | | | fwd | 1.4e-15 | 0.059 | 0.335 | 0.053 | 0.393 |
| variable_element | haunched `VariableElement` model | 9 | FD | rev | 1.6e-09 | 0.161 | — | 0.202 | — |
| | | | | fwd | 1.6e-09 | 0.274 | — | 0.319 | — |
| diaphragm | frame with a rigid diaphragm | 2 | Zygote | rev | 7.6e-15 | 0.149 | 0.447 | 0.134 | 0.484 |
| | | | | fwd | 6.4e-15 | 0.086 | 0.428 | 0.092 | 0.465 |
| shells | shell thickness (`ThicknessVariable`) | 2 | Zygote | rev | 0.0e+00 | 6.68 | 7.74 | 6.66 | 9.52 |
| | | | | fwd | 1.0e-14 | 1.32 | 9.52 | 1.1 | 8.5 |
| shells_tri | shell thickness on the same wall, triangulated (ANDES-OPT + MITC3+) | 2 | Zygote | rev | 1.9e-16 | 5.88 | 7.97 | 5.09 | 9.22 |
| | | | | fwd | 3.1e-15 | 1.14 | 8.25 | 0.98 | 9.51 |

First-call (compile) time per process: 27–60 s for most cases, up to 82 s for shells in reverse mode.

## Reading the table

- **Reverse mode is the default for gradients.** It is faster than Zygote in every case, by 1.2× (shells) to 5.4× (frame line loads on 1.12). The 1.3 Zygote rules (`AsapFrameADExt`, and the linear-cost group rules in `AsapChainRulesExt`) narrowed the gap; on these small models the group rules matter little, their gain grows with model size.
- **Forward mode costs one pass per variable**, so it wins only for a handful of variables
  (shells: 1.3 ms against 6.7 ms in reverse mode) and loses badly at 512 (285–538 ms). For
  Jacobians of many outputs with respect to few inputs, forward mode is the right tool.
- **Zygote cannot differentiate `VariableElement` models** ("Mutating arrays is not supported").
  Enzyme can. This predates 1.3.
- **Shells** (quadrilaterals and, since 1.4, triangles) run through the bridged rules in
  `AsapShellEnzymeExt` in reverse mode, and natively in forward mode.

## What changed since the August 2026 status

In August (Enzyme 0.13.199), forward mode aborted on Julia 1.12 with an `AdjointGenerator.h`
assertion, threw `ReadOnlyMemoryError` at spaceframe scale on 1.11, and an older Enzyme hung in
compilation. With 0.13.221, forward mode works on 1.12, and spaceframe-scale forward mode works on
1.12 and 1.13. Julia 1.11 was not rerun, so its two failures are untested.

These results are from the final 1.3 code (rerun 2026-11-24, after the linear-cost Zygote rules and
the split shell rules). The matrix has caught two regressions that the Asap suite missed, both
reverted before commit:

- a pure-path refactor (gather + `map` over abstractly typed arrays): Enzyme reverse threw, and
  forward returned wrong derivatives silently;
- factoring the frame stiffness lines of `_group_entries` into a helper function: Enzyme reverse
  failed type analysis on `frame_joints` ("Cannot deduce type of load").

Rerun the matrix after any change to code on the pure path, however small. The shell rules
bridged into Enzyme (`_shell_K_plain`, `_shell_fef_plain`) got faster in the same change, so
shell reverse mode went from 8.5 to 6.4 ms.

One Asap fix was needed for 1.12. `AsapEnzymeExt` and `AsapShellEnzymeExt` register their rules with
`@eval` in `__init__`. While `AsapShellEnzymeExt` precompiled, it loaded `AsapEnzymeExt`, whose
`__init__` then evaluated into a closed module, a hard error on 1.12. Both `__init__`s now return
early while output is being generated (`ccall(:jl_generating_output, Cint, ()) == 1`). The rules
are registered when the extensions are actually loaded.

Keep wrapping Enzyme runs in an external `timeout`; its failure modes have included silent
compile-time hangs.
