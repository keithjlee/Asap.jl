# Shell element derivation dossier (S0)

> Formulations in Asap notation, transcribed from the primary sources. Every equation that came
> from a paper carries its source equation number so it can be re-checked.
> **Status:** membrane (§2) transcribed from the ITW report (verified page by page,
> 2026-11-23). Plate (§3) is written from the standard MITC4 formulation. Its tying-point
> scheme should be re-checked against Bathe & Dvorkin (1985) when that PDF is pulled
> (open item below).

## 1. Conventions

**Natural coordinates and node order (Asap).** A quad has four corner nodes, ordered
counter-clockwise about the element normal:

| Asap node a | (ξₐ, ηₐ) | ITW node (Fig. 1) |
|---|---|---|
| 1 | (−1, −1) | 3 |
| 2 | (+1, −1) | 4 |
| 3 | (+1, +1) | 1 |
| 4 | (−1, +1) | 2 |

Both orderings are counter-clockwise; ours is ITW's shifted by two. Every ITW formula below is
written per **edge**, which makes it invariant under a cyclic relabelling.

Edge k joins nodes k → k+1 (with node 5 ≡ node 1), and its midside serendipity function is:

| Edge k | Nodes | Midside point | NSₖ(ξ, η) |
|---|---|---|---|
| 1 | 1 → 2 | (0, −1) | ½(1 − ξ²)(1 − η) |
| 2 | 2 → 3 | (+1, 0) | ½(1 + ξ)(1 − η²) |
| 3 | 3 → 4 | (0, +1) | ½(1 − ξ²)(1 + η) |
| 4 | 4 → 1 | (−1, 0) | ½(1 − ξ)(1 − η²) |

These are ITW (3.7)–(3.8). Corner functions are Nₐ = ¼(1 + ξₐξ)(1 + ηₐη) (ITW 3.3). The bubble
is NB = (1 − ξ²)(1 − η²) (ITW 3.9).

**Local frame.** The four nodes give global positions Xₐ. Then:

- e₃ = normalize((X₃ − X₁) × (X₄ − X₂)), the normal from the diagonals, which is well defined
  for warped quads;
- e₁ = normalize(projection of (X₂ − X₁) onto the plane ⊥ e₃), so edge 1 sets the local x-axis
  (an optional material-axis vector may replace it later);
- e₂ = e₃ × e₁;
- R = [e₁ e₂ e₃]ᵀ (rows = local axes), so v_local = R·v_global;
- local 2D nodal coordinates: (xₐ, yₐ) = first two components of R·(Xₐ − X̄), with X̄ the
  centroid. Warped nodes are projected onto the mean plane. The warping ratio
  |R·(Xₐ − X̄)|₃ / √area is checked, with a warning above a threshold (MacNeal–Harder style).
  A warping correction is deferred (S7).

**Nodal DOFs.** Asap gives every node 6 global slots: (uₓ, u_y, u_z, θₓ, θ_y, θ_z). Per node,
the local DOFs are (u, v, w, θx, θy, ψ) = (R·u_global, R·θ_global). Here θx, θy are rotations
about the local x and y axes, and ψ ≡ θz is the **drilling rotation**, counter-clockwise about
e₃. The global element matrix is K_g = Tᵀ K_l T, with T = blockdiag(R, R, …) over the 8 3×3
blocks. That is applied blockwise, never as a 24×24 product.

**Split.** The membrane uses (u, v, ψ) and the plate uses (w, θx, θy). There is no
membrane–plate coupling in a flat isotropic element (or any symmetric section). An
eccentric/layered section would add a B coupling block, but that is out of scope.

**Section rigidities.** Material E, ν, thickness t, μ = G = E/(2(1 + ν)):

| Symbol | Meaning | Value | Dimension |
|---|---|---|---|
| Dm | membrane | t · E/(1 − ν²) · [1 ν 0; ν 1 0; 0 0 (1 − ν)/2] | [force/length] |
| Db | bending | t³/12 · E/(1 − ν²) · [same matrix] | [force·length] |
| Ds | transverse shear | κ·G·t·I₂, κ = 5/6 | [force/length] |
| γₜ | drilling penalty | γ·t, with γ = μ (ITW §2; Hughes et al. 1989), so γₜ = Dm[3,3] | [force/length] |

ITW Table 5 shows the solution is insensitive to γ/μ ∈ [10⁻³, 10³].

## 2. Membrane: Ibrahimbegović–Taylor–Wilson drilling quad, D-type

Source: Ibrahimbegović, Taylor & Wilson, UCB/SEMM-89/17 (1989) = IJNME 30:445–457 (1990).
We implement the **displacement-type** ("D-type") formulation, ITW Problem (Dʰ), eq. (3.22).

**Functional** (ITW 2.14), per unit thickness and then multiplied by t:

  Π(u, ψ) = ½∫ sym∇u : C : sym∇u dΩ + ½ γ ∫ (skew∇u − ψ)² dΩ − ∫ u·f dΩ

with skew∇u ≡ ½(∂v/∂x − ∂u/∂y), the scalar infinitesimal rotation.

**Rotation field** (ITW 3.2): ψʰ = Σₐ Nₐ ψₐ (bilinear).

**Displacement field** (ITW 3.4, Allman-type, rewritten per edge):

  uʰ = Σₐ Nₐ uₐ + Σₖ (lₖ/8)(ψ_{k+1} − ψ_k) nₖ NSₖ + NB·Δu_b

where, for edge k from node k to node k+1:

- lₖ = |x_{k+1} − x_k|;
- nₖ is the outward unit normal (ITW 3.5). For counter-clockwise nodes,
  lₖ nₖ = (y_{k+1} − y_k, −(x_{k+1} − x_k)), so the products lₖnₖ are polynomial in the
  coordinates;
- the bubble amplitude Δu_b ∈ ℝ² is condensed statically (Wilson 1974).

Equivalently, node a's drilling shape vector is

  dₐ(ξ, η) = (1/8)[ l_{a−1} n_{a−1} NS_{a−1} − lₐ nₐ NSₐ ]

where edge a−1 ends at node a (coefficient +ψₐ) and edge a starts at node a (coefficient −ψₐ).
This is ITW (3.13)–(3.14) with the J, K, L, M index gymnastics removed. Then
u_drill = Σₐ dₐ ψₐ.

**Discrete operators** (ITW 3.11–3.17), per node a:

- symmetric strain: ε = Bₐ uₐ + Gₐ ψₐ, with Bₐ the standard [Nₐ,ₓ 0; 0 Nₐ,ᵧ; Nₐ,ᵧ Nₐ,ₓ]
  and Gₐ = [dₐ,ₓ⁽¹⁾; dₐ,ᵧ⁽²⁾; dₐ,ᵧ⁽¹⁾ + dₐ,ₓ⁽²⁾] (ITW 3.12–3.13);
- skew minus rotation: skew∇u − ψ = bₐ·uₐ + gₐ ψₐ, with
  bₐ = (−½ Nₐ,ᵧ, ½ Nₐ,ₓ) (ITW 3.16) and gₐ = ½(dₐ,ₓ⁽²⁾ − dₐ,ᵧ⁽¹⁾) − Nₐ (ITW 3.17);
- bubble strain columns Bb = [NB,ₓ 0; 0 NB,ᵧ; NB,ᵧ NB,ₓ]. The bubble's skew contribution
  vanishes at the element centre (∇NB(0, 0) = 0), so it adds nothing to P under 1-point
  quadrature (ITW p. 11).

**Element matrices:**

- K = ∫ [B G Bb]ᵀ Dm [B G Bb] dΩ with **3×3 Gauss** (ITW p. 10), size 14×14 over
  (u, v, ψ)×4 plus the bubble's 2 DOFs. Condense the bubble:
  K̂ = K_aa − K_ab K_bb⁻¹ K_ba.
- P = γₜ ∫ [b g]ᵀ[b g] dΩ with **1-point Gauss** (ITW 3.23, p. 10).
- K_membrane = K̂ + P (ITW 3.24), size 12×12.

Full integration of K combined with P removes the spurious zero-energy modes, so no other
device is needed (ITW p. 10). Allman / Taylor–Simo elements have a spurious constant-ψ mode;
this one does not (ITW §4.1).

**Verification targets** (ITW §4):

- **Patch test** (Fig. 2): skewed single element under uniform tension, minimal restraints, must
  give exact constant stress.
- **Simple beam** (Fig. 3, Table 1): E = 100, ν = 0.3, l = 10, h = 1, one row of 6 elements,
  minimal restraints, no drilling restraint.

  | Load case | Mesh | Vertical displacement | End rotation |
  |---|---|---|---|
  | LC1 (end couple, nodal forces P = 1) | regular | 1.5 | 0.6 |
  | LC1 | distorted | 1.14045 | 0.57247 |
  | LC2 (M = 0.5 at the drilling DOFs at both ends) | regular | 1.5 | 0.62070 |
  | LC2 | distorted | 1.39200 | 0.49508 |

  All values are D-type. The beam-theory exact solution is 1.5 / 0.6. **LC2 directly tests
  the drilling sign convention.**
- **Short cantilever** (Fig. 4, Table 2): E = 30000, ν = 0.25, P = 40, l = 48, h = 12 (unit
  thickness), parabolic shear at the tip. Exact tip deflection
  Pl³/(3EI) + (4 + 5ν)Pl/(2Eh) = 0.3553. Fixed end: one pin and two rollers (horizontal-only
  restraint).

  | Mesh | D-type |
  |---|---|
  | 4×1 | 0.3445 |
  | 8×2 | 0.3504 |
  | 16×4 | 0.3543 |
  | 4×1 distorted (MacNeal–Harder) | 0.3065 |

- **Cook's membrane** (Fig. 5, Table 3): E = 1, ν = 1/3, P = 1 distributed on the right edge,
  l = 48, h1 = 16, h2 = 44 (left edge 44 tall, right edge from y = 44 to y = 60).
  Reference tip deflection 23.91.

  | Mesh | D-type |
  |---|---|
  | 1×1 | 14.065 |
  | 2×2 | 20.682 |
  | 4×4 | 22.984 |
  | 8×8 | 23.626 |

- **Pinched hemisphere with 18° hole** (Fig. 6, Table 4): ITW pair the membrane with **DKQ**,
  so these values are a sanity range only for our MITC4 pairing. E = 6.825e7, ν = 0.3,
  R = 10, t = 0.04, P = 1 per quadrant. Reference 0.094 (MacNeal–Harder) / 0.093 (Simo et al.).

  | Mesh | ITW D-type |
  |---|---|
  | 4×4 | 0.0875 |
  | 8×8 | 0.0937 |
  | 16×16 | 0.0935 |

### 2.1 Implementation findings (S2, 2026-11-23) and the enhanced variant

**Verified properties of the implementation.**

- Exactly 3 zero eigenvalues (u, v and in-plane rotation with ψ = θ). A constant ψ is not a
  zero mode.
- Both variants pass the distorted-patch test, including ψ = ω, to 1e-14.
- The strain operator matches finite differences of the interpolated field.
- Displacements are continuous across shared edges (checked numerically).
- ½aᵀKa equals an independent 12×12-Gauss integration of the field energy: exact on
  rectangles, 3e-5 on a distorted quad from 3×3 Gauss on a rational integrand.

**Why the ITW tables are not reproduced digit for digit.**

1. **Drilling DOFs at supports.** On an Allman element, u = v = 0 at the nodes of an edge does
   *not* clamp the edge. With ψ free, the edge term (l/8)(ψ_K − ψ_J)n lets it bulge
   quadratically between the nodes.
   - So "clamped" must also fix ψ, which Asap's `:fixed` does.
   - With ψ free, the clamped cantilever overshoots the converged value by up to 7% (8×2) and
     converges non-monotonically.
   - ITW's figures show point supports, and the treatment of ψ there is not stated.
2. **Cook's reference point.** The literature value 23.96 (ITW: 23.91; Bergan & Fellipa 1985:
   23.9) is the vertical displacement at the **right-edge midpoint (48, 52)**.
   - The top corner converges to about 25.1 (conforming Q4 and ITW agree).
   - ITW's own 1×1 entry (14.065) can't be a midpoint value, so their measurement point and
     boundary treatment are ambiguous.
3. **Pure bending on a regular mesh (ITW Table 1).** ITW report exactly 1.5 / 0.6. The
   interpolation as transcribed gives 1.4615 / 0.5801 (ν = 0.3).
   - Allman edge modes span only the *normal* edge quadratics, and one bubble cannot produce
     the Poisson strain ε_yy = νκy uniformly.
   - Without the bubble the result is exactly 1.5·(1 − ν²) = 1.365, the classic Poisson
     locking. This is an analytical property of the transcribed fields, not a coding error.
     ITW's code presumably had more than the paper describes.

**Enhanced variant (the default; `enhanced = false` gives the element as published).** Add the
two *complementary* tangential incompatible modes, u ∝ 1 − ξ² and v ∝ 1 − η².

- Derivatives use the Taylor–Beresford–Wilson (1976) correction: centre Jacobian, scaled by
  detJ₀/detJ, so the distorted patch test still passes.
- Both modes are condensed together with the bubble. They carry no rotation into the 1-point
  penalty.
- Adding all four Wilson modes instead (QM6) duplicates Allman's normal quadratics and makes ψ
  ill-determined (end rotation 23 instead of 0.6), so only the complementary pair is used.

**Benchmark comparison** (fully clamped; ratio to a converged reference cross-checked against a
256×64 conforming-Q4 run, 0.35601):

| Case | Published ITW | Enhanced |
|---|---|---|
| Pure bending, regular mesh | 1.4615 / 0.5801 | **1.5 / 0.6 (exact)** |
| Cantilever l/h = 4, 4×1 · 8×2 · 16×4 | 0.955 · 0.987 · 0.995 | **0.981 · 1.000 · 0.998** |
| Cantilever, distorted 4×1 | 0.939 | **0.963** |
| Squat wall l/h = 1, 1×1 · 2×2 · 4×4 | 0.757 · 0.926 · 0.973 | **0.790 · 0.969 · 0.983** |
| Cook midpoint 2 · 4 · 8 (→ 23.96) | 21.62 · 23.17 · 23.66 | **22.03 · 23.32 · 23.70** |

For scale, plain Q4 gives 0.68 on the 4×1 cantilever. A wall meshed two elements wide is now
within 3% in-plane, which was the point of choosing a drilling membrane.

## 3. Plate: MITC4 (Dvorkin–Bathe)

Sources: Dvorkin & Bathe (1984), Engineering Computations 1:77–88; Bathe & Dvorkin (1985),
IJNME 21:367–383.

**Kinematics** (Reissner–Mindlin, linear): the normal rotates by β = (βx, βy), with in-plane
displacement through the thickness z·β.

- With right-hand nodal rotations about local axes: βx = θy and βy = −θx.
- w, θx, θy are interpolated bilinearly with Nₐ.
- Curvatures: κ = [βx,ₓ ; βy,ᵧ ; βx,ᵧ + βy,ₓ]. Bending energy ½∫ κᵀ Db κ dΩ, **2×2 Gauss**.
- Transverse shear, Cartesian: γ = (w,ₓ + βx, w,ᵧ + βy).

**MITC assumed transverse shear.** Use covariant components in natural coordinates:

  e_ξ = w,ξ + β·x,ξ   and   e_η = w,η + β·x,η,  with x,ξ = (x,ξ, y,ξ)

evaluated with the displacement interpolation at the tying points:

- e_ξ at A = (0, −1) and C = (0, +1), the midpoints of edges 1 and 3;
- e_η at B = (+1, 0) and D = (−1, 0), the midpoints of edges 2 and 4.

These are interpolated as

  ẽ_ξ(ξ, η) = ½(1 − η) e_ξᴬ + ½(1 + η) e_ξᶜ
  ẽ_η(ξ, η) = ½(1 + ξ) e_ηᴮ + ½(1 − ξ) e_ηᴰ

and mapped to Cartesian components with γ = J⁻¹·(ẽ_ξ, ẽ_η), where J = [x,ξ y,ξ; x,η y,η] so
that (∂/∂ξ, ∂/∂η) = J·(∂/∂x, ∂/∂y). Shear energy ½∫ γᵀ Ds γ dΩ, **2×2 Gauss**.

- Each covariant component is linear along its edge, so the tying values are exact for any
  bilinear geometry. This is what removes shear locking: in the thin limit the Kirchhoff
  constraint is enforced only at the tying points.
- MITC4 has no spurious zero-energy modes with full 2×2 integration, and passes the patch test.

**Implemented in S3 (2026-11-23), results:**
- exactly 3 rigid modes for the plate and 6 for the full shell in any orientation;
- the constant-curvature patch test passes on a distorted patch for thin (t = 0.001) and thick
  (0.05) plates;
- square plate, simply supported, uniform load:
  0.0039728 · 0.0040466 · 0.0040653 (4×4 · 8×8 · 16×16, reference 0.00406);
- clamped plate: 0.0012641 at 16×16 (reference 0.00126);
- centre point load: 0.011612 at 16×16 (reference 0.0116);
- thickness sweep (8×8, a/t = 10² · 10³ · 10⁴): 0.0040466 · 0.0040415 · 0.0040414, so no
  locking; a/t = 10 gives 0.0044731 (Mindlin shear flexibility);
- Scordelis–Lo, lumped gravity: 0.3174 · 0.3042 · 0.3017 · 0.3014 (4 · 8 · 16 · 32 mesh;
  reference 0.3024, and flat facets converge to about 0.301).

**Still open:** pull Bathe & Dvorkin (1985) and confirm the tying-point
interpolation and sign conventions verbatim, and record their benchmark values. That is the
"derivations verifiable against primary literature" standard.

**Verification targets (S3)** (references to be pulled in S0/S3):

- Rigid-body modes: the full 24×24 element has exactly 6 zero eigenvalues (any geometry),
  including warped.
- Constant-curvature plate patch test (MacNeal–Harder patch geometry).
- Simply supported square plate under uniform load q, Kirchhoff limit: w_c = 0.00406 qa⁴/D
  (ν = 0.3). Clamped: w_c = 0.00126 qa⁴/D. Central point load P, simply supported:
  w_c = 0.0116 Pa²/D (Timoshenko & Woinowsky-Krieger).
- Thin-to-thick sweep (a/t = 10 … 10⁴) to confirm there is no locking.
- MacNeal–Harder (1985):
  - Scordelis–Lo roof: 0.3024;
  - pinched cylinder with diaphragms: 1.8248e-5;
  - twisted beam: 5.424e-3 in-plane, 1.754e-3 out-of-plane;
  - pinched hemisphere: 0.094;
  - straight cantilever with rectangular, trapezoidal and parallelogram meshes.

## 4. Stress resultants (S4)

At any (ξ, η), from the local element displacements:

- N = Dm · ε (membrane, including the condensed bubble recovered from Δu_b = −K_bb⁻¹ K_ba a);
- M = Db · κ;
- Q = Ds · γ̃ (MITC field);
- local nodal forces f_l = K_l a_l + q_l (with q_l the consistent element-load vector), used
  for section cuts (§ plan 4). They are exact for equilibrium: Σ over a free body of nodal
  forces = applied loads.

## 5. Rigid diaphragm constraint (S5b)

For a diaphragm with unit normal n (default global Z) and master node m:

- each slave s has in-plane translations and the rotation about n tied to the master:
  u_s∥ = u_m∥ + θ_m,n · n × (x_s − x_m) and θ_s,n = θ_m,n;
- the out-of-plane translation and the two out-of-plane rotations stay independent;
- implemented as a transformation u = T û, applied **per element** (K̂ₑ = Tₑᵀ Kₑ Tₑ, where Tₑ is
  a pure function of positions, so it is AD-transparent and the frozen-pattern / nzmap
  machinery still holds);
- loads P̂ = Tᵀ P, and recovery u = T û;
- the master is an existing node or an auto-generated node at the diaphragm centroid.

Details to be written up in S5b before implementing.
