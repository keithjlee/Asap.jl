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
| truss_area | member areas, planar truss | 14 | Zygote | rev | 8.8e-16 | 0.083 | 0.286 | 0.09 | 0.414 |
| | | | | fwd | 2.8e-15 | 0.119 | 0.278 | 0.114 | 0.393 |
| truss_shape | top-chord node heights | 3 | Zygote | rev | 6.4e-16 | 0.09 | 0.311 | 0.1 | 0.429 |
| | | | | fwd | 1.1e-15 | 0.059 | 0.35 | 0.054 | 0.399 |
| spaceframe | member areas, 8 × 8 space truss | 512 | Zygote | rev | 1.0e-14 | 0.861 | 1.18 | 0.877 | 1.25 |
| | | | | fwd | 1.5e-14 | 168 | 1.92 | 584 | 35.9 |
| frame_lineload | frame section with line load + self-weight (design-dependent FEFs) | 3 | Zygote | rev | 7.7e-15 | 0.119 | 0.696 | 0.125 | 0.968 |
| | | | | fwd | 4.1e-15 | 0.109 | 0.651 | 0.112 | 1.11 |
| frame_joints | semi-rigid connection stiffness (`JointVariable`) | 1 | Zygote | rev | 1.9e-16 | 0.512 | 0.973 | 0.547 | 1.31 |
| | | | | fwd | 3.8e-16 | 0.468 | 0.967 | 0.452 | 1.17 |
| springs | spring-support stiffness | 1 | Zygote | rev | 8.1e-16 | 0.091 | 0.371 | 0.112 | 0.521 |
| | | | | fwd | 1.4e-15 | 0.065 | 0.397 | 0.063 | 0.576 |
| variable_element | haunched `VariableElement` model | 9 | FD | rev | 1.6e-09 | 0.126 | — | 0.152 | — |
| | | | | fwd | 1.6e-09 | 0.239 | — | 0.281 | — |
| diaphragm | frame with a rigid diaphragm | 2 | Zygote | rev | 7.6e-15 | 0.244 | 0.545 | 0.14 | 0.738 |
| | | | | fwd | 6.4e-15 | 0.079 | 0.508 | 0.091 | 0.763 |
| shells | shell thickness (`ThicknessVariable`) | 2 | Zygote | rev | 0.0e+00 | 8.47 | 13.5 | 8.7 | 15.3 |
| | | | | fwd | 1.0e-14 | 1.01 | 13.3 | 1.1 | 14.4 |

First-call (compile) time per process: 26–54 s for most cases, up to 95 s for shells in reverse mode.

## Reading the table

- **Reverse mode is the default for gradients.** It is faster than Zygote in every case, by 1.4× (spaceframe) to 7.7× (frame line loads on 1.12). Zygote frame gradients use the 1.3 frame rules (`AsapFrameADExt`), which narrowed the gap.
- **Forward mode costs one pass per variable**, so it wins only for a handful of variables
  (shells: 1.1 ms against 8 ms in reverse mode) and loses badly at 512 (168–584 ms). For
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

These results are from the final 1.3 code (rerun 2026-11-24). An intermediate pure-path refactor
passed the Asap suite but failed this matrix (Enzyme reverse threw; forward returned wrong
derivatives silently) and was reverted: rerun the matrix after any pure-path change.

One Asap fix was needed for 1.12. `AsapEnzymeExt` and `AsapShellEnzymeExt` register their rules with
`@eval` in `__init__`. While `AsapShellEnzymeExt` precompiled, it loaded `AsapEnzymeExt`, whose
`__init__` then evaluated into a closed module, a hard error on 1.12. Both `__init__`s now return
early while output is being generated (`ccall(:jl_generating_output, Cint, ()) == 1`). The rules
are registered when the extensions are actually loaded.

Keep wrapping Enzyme runs in an external `timeout`; its failure modes have included silent
compile-time hangs.
