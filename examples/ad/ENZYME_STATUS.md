# Enzyme across Asap — status (2026-11-24)

**Enzyme works across Asap, in reverse and forward mode, on Julia 1.13 and 1.12.** All 18 cases
(9 scenario classes × 2 modes) pass on both Julia versions. Gradients agree with Zygote to
1e-14 or better. On `VariableElement` models, which Zygote cannot differentiate, they agree with
central finite differences to 1.6e-9, the finite-difference truncation error.

Versions: Enzyme 0.13.221, Zygote 0.7.20, DifferentiationInterface 0.7.22; Asap and AsapOptim on
branch `shells`. Settings as everywhere in Asap: `set_runtime_activity(Reverse / Forward)` with a
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
| truss_area | member areas, planar truss | 14 | Zygote | rev | 8.8e-16 | 0.094 | 0.301 | 0.093 | 0.519 |
| | | | | fwd | 2.8e-15 | 0.122 | 0.317 | 0.126 | 0.477 |
| truss_shape | top-chord node heights | 3 | Zygote | rev | 6.4e-16 | 0.091 | 0.301 | 0.101 | 0.418 |
| | | | | fwd | 1.1e-15 | 0.054 | 0.328 | 0.056 | 0.462 |
| spaceframe | member areas, 8 × 8 space truss | 512 | Zygote | rev | 1.0e-14 | 0.909 | 1.22 | 1.11 | 1.4 |
| | | | | fwd | 1.5e-14 | 200 | 2.48 | 546 | 1.22 |
| frame_lineload | frame section with line load + self-weight (design-dependent FEFs) | 3 | Zygote | rev | 7.7e-15 | 0.122 | 0.521 | 0.312 | 0.904 |
| | | | | fwd | 4.1e-15 | 0.094 | 0.56 | 0.113 | 0.776 |
| frame_joints | semi-rigid connection stiffness (`JointVariable`) | 1 | Zygote | rev | 1.9e-16 | 0.484 | 0.828 | 0.602 | 1.17 |
| | | | | fwd | 3.8e-16 | 0.451 | 0.812 | 0.541 | 1.16 |
| springs | spring-support stiffness | 1 | Zygote | rev | 2.0e-16 | 0.09 | 0.278 | 0.104 | 0.355 |
| | | | | fwd | 1.4e-15 | 0.054 | 0.311 | 0.057 | 0.396 |
| variable_element | haunched `VariableElement` model | 9 | FD | rev | 1.6e-09 | 0.122 | — | 0.156 | — |
| | | | | fwd | 1.6e-09 | 0.264 | — | 0.456 | — |
| diaphragm | frame with a rigid diaphragm | 2 | Zygote | rev | 7.6e-15 | 0.128 | 0.386 | 0.127 | 0.456 |
| | | | | fwd | 6.4e-15 | 0.082 | 0.389 | 0.14 | 0.472 |
| shells | shell thickness (`ThicknessVariable`) | 2 | Zygote | rev | 0.0e+00 | 6.38 | 8.08 | 6.92 | 8.95 |
| | | | | fwd | 1.0e-14 | 1.11 | 8.07 | 1.2 | 9.09 |

First-call (compile) time per process: 28–59 s for most cases, up to 84 s for shells in reverse mode.

## Reading the table

- **Reverse mode is the default for gradients.** It is faster than Zygote in every case, by 1.3× (shells, spaceframe) to 5.6× (truss areas on 1.12). The 1.3 Zygote rules (`AsapFrameADExt`, and the linear-cost group rules in `AsapChainRulesExt`) narrowed the gap; on these small models the group rules matter little, their gain grows with model size.
- **Forward mode costs one pass per variable**, so it wins only for a handful of variables
  (shells: 1.1 ms against 6.4 ms in reverse mode) and loses badly at 512 (200–546 ms). For
  Jacobians of many outputs with respect to few inputs, forward mode is the right tool.
- **Zygote cannot differentiate `VariableElement` models** ("Mutating arrays is not supported").
  Enzyme can. This predates 1.3.
- **Shells** run through the bridged rules in `AsapShellEnzymeExt` in reverse mode, and natively
  in forward mode.

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
