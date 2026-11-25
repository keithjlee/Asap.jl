#=
Flat 3-node shell kernels (triangles) — pure functions of the local node
coordinates and the section rigidities, like the quad kernels in shell.jl.

  membrane: ANDES-OPT drilling triangle (Felippa 2003, CMAME 192:2125–2168),
            the optimal instance of the LST-3/9R ANDES template
  plate:    MITC3+ (Lee, Lee & Bathe 2014, Comput. Struct. 138:12–23),
            Reissner–Mindlin with a condensed cubic rotation bubble and the
            MITC3+ assumed transverse shear field

Node order: counter-clockwise about the element normal. Local per-node DOFs
as for the quad: membrane (u, v, ψ), plate (w, θx, θy), with ψ = θz the
drilling rotation and the plate using the normal's rotations
β = (βx, βy) = (θy, −θx).
=#

# ─────────────────────────────────────────────────────────────────────────────
# Membrane: ANDES-OPT
# ─────────────────────────────────────────────────────────────────────────────

# OPT signature (Felippa 2003, eq 42): α_b and β₁…β₉; β₀ follows the material
const _ANDES_ALPHA_B = 1.5
const _ANDES_BETA = (1.0, 2.0, 1.0, 0.0, 1.0, -1.0, -1.0, -1.0, -2.0)

"""
    andes_beta0(A) -> T

Higher-order stiffness scaling `β₀` of the ANDES-OPT triangle for the membrane
rigidity `A` (Felippa 2003, eq 44): `β₀ = 2/Ῡ − 3/2`, with `Ῡ` the average over
all directions φ of `E₁₁(φ)·C₁₁(φ)` (rotated stiffness × compliance), clamped
below at 0.01 to keep the element rank sufficient. Isotropic: `½(1 − 4ν²)`.

`E₁₁(φ)·C₁₁(φ)` is a trigonometric polynomial of degree 8 in φ, so its
average is exact with 16 equally spaced samples.
"""
function andes_beta0(A::AbstractMatrix)
    C = inv(SMatrix{3,3}(A))
    Υ = sum(1:16) do k
        φ = (k - 1) * (π / 8)
        c, s = cos(φ), sin(φ)
        E11 = A[1, 1] * c^4 + 4A[1, 3] * c^3 * s + (2A[1, 2] + 4A[3, 3]) * c^2 * s^2 +
              4A[2, 3] * c * s^3 + A[2, 2] * s^4
        C11 = C[1, 1] * c^4 + 2C[1, 3] * c^3 * s + (2C[1, 2] + C[3, 3]) * c^2 * s^2 +
              2C[2, 3] * c * s^3 + C[2, 2] * s^4
        E11 * C11
    end / 16
    return max(2 / Υ - 3 / 2, oftype(Υ, 1 / 100))
end

"""
    andes_membrane_stiffness(xy, A) -> SMatrix{9,9}

Local ANDES-OPT membrane stiffness of a flat triangle over the DOFs
`(u, v, ψ)` of nodes 1…3 (Felippa 2003): `K = K_b + K_h`, with the basic
stiffness `K_b = LALᵀ/area` (force lumping matrix `L`, α_b = 3/2) and the
higher-order stiffness `K_h = T̃ᵀ K_θ T̃` built from deviatoric natural strains
along the sides (β₁…β₉ of the optimal signature, β₀ from [`andes_beta0`](@ref)).

- `xy`: local in-plane node coordinates, counter-clockwise [length]
- `A`: membrane rigidity (thickness-integrated, `N = A·ε`) [force/length]

Exact for in-plane pure bending of rectangular two-triangle mesh units of any
aspect ratio; passes the individual element test (patch test) for any β.
"""
function andes_membrane_stiffness(xy, A::AbstractMatrix)
    E = SMatrix{3,3}(A)
    L, Ar, Tθu, Te, Q1, Q2, Q3 = _andes_parts(xy)
    Kb = L * E * L' / Ar
    Enat = Te' * E * Te
    Q4, Q5, Q6 = (Q1 + Q2) / 2, (Q2 + Q3) / 2, (Q3 + Q1) / 2
    Kθ = (3 / 4) * andes_beta0(E) * Ar * (Q4' * Enat * Q4 + Q5' * Enat * Q5 + Q6' * Enat * Q6)
    return Kb + Tθu' * Kθ * Tθu
end

# geometry-only parts of the ANDES template (Felippa 2003 module, Fig. 10):
# force-lumping matrix L (9×3, h = 1, eq 22), signed area, hierarchical-rotation
# map Tθu (eq 16), natural-to-Cartesian strain map Te (eq 12) and Q₁…Q₃ (eq 29)
function _andes_parts(xy)
    (x1, y1), (x2, y2), (x3, y3) = xy[1], xy[2], xy[3]
    x12, x23, x31 = x1 - x2, x2 - x3, x3 - x1
    y12, y23, y31 = y1 - y2, y2 - y3, y3 - y1
    x21, x32, x13 = -x12, -x23, -x31
    y21, y32, y13 = -y12, -y23, -y31
    Ar = (y21 * x13 - x21 * y13) / 2                  # signed area (> 0 counter-clockwise)
    A2, A4 = 2Ar, 4Ar
    ab = _ANDES_ALPHA_B
    L = SMatrix{9,3}(
        y23, 0, ab * y23 * (y13 - y21) / 6, y31, 0, ab * y31 * (y21 - y32) / 6, y12, 0, ab * y12 * (y32 - y13) / 6,
        0, x32, ab * x32 * (x31 - x12) / 6, 0, x13, ab * x13 * (x12 - x23) / 6, 0, x21, ab * x21 * (x23 - x31) / 6,
        x32, y23, ab * (x31 * y13 - x12 * y21) / 3, x13, y31, ab * (x12 * y21 - x23 * y32) / 3,
        x21, y12, ab * (x23 * y32 - x31 * y13) / 3) / 2
    Tθu = SMatrix{3,9}(
        x32, x32, x32, y32, y32, y32, A4, 0, 0,
        x13, x13, x13, y13, y13, y13, 0, A4, 0,
        x21, x21, x21, y21, y21, y21, 0, 0, A4) / A4
    LL21, LL32, LL13 = x21^2 + y21^2, x32^2 + y32^2, x13^2 + y13^2
    Te = SMatrix{3,3}(
        y23 * y13 * LL21, x23 * x13 * LL21, (y23 * x31 + x32 * y13) * LL21,
        y31 * y21 * LL32, x31 * x21 * LL32, (y31 * x12 + x13 * y21) * LL32,
        y12 * y32 * LL13, x12 * x32 * LL13, (y12 * x23 + x21 * y32) * LL13) / (Ar * A4)
    b1, b2, b3, b4, b5, b6, b7, b8, b9 = _ANDES_BETA
    c = A2 / 3
    Q1 = SMatrix{3,3}(b1 / LL21, b4 / LL32, b7 / LL13, b2 / LL21, b5 / LL32, b8 / LL13, b3 / LL21, b6 / LL32, b9 / LL13) * c
    Q2 = SMatrix{3,3}(b9 / LL21, b3 / LL32, b6 / LL13, b7 / LL21, b1 / LL32, b4 / LL13, b8 / LL21, b2 / LL32, b5 / LL13) * c
    Q3 = SMatrix{3,3}(b5 / LL21, b8 / LL32, b2 / LL13, b6 / LL21, b9 / LL32, b3 / LL13, b4 / LL21, b7 / LL32, b1 / LL13) * c
    return L, Ar, Tθu, Te, Q1, Q2, Q3
end

# Cartesian membrane strain (εxx, εyy, γxy) at triangular coordinates
# (ζ₁, ζ₂, ζ₃) = (1−r−s, r, s) from the local membrane DOFs a = (u, v, ψ)×3:
# e = Lᵀa/area + Te·β₀ᵉ·(Q₁ζ₁ + Q₂ζ₂ + Q₃ζ₃)·Tθu·a (Felippa eq 49), with the
# recommended recovery scaling β₀ᵉ = 3/2. At the centroid the higher-order term
# vanishes (the OPT element is energy-orthogonal: Q₁ + Q₂ + Q₃ = 0).
function _andes_strain(xy, a::SVector{9}, r, s)
    L, Ar, Tθu, Te, Q1, Q2, Q3 = _andes_parts(xy)
    Qζ = (1 - r - s) * Q1 + r * Q2 + s * Q3
    return L' * a / Ar + Te * ((3 / 2) * (Qζ * (Tθu * a)))
end

# ─────────────────────────────────────────────────────────────────────────────
# Plate: MITC3+
# ─────────────────────────────────────────────────────────────────────────────

# 7-point degree-5 rule on the reference triangle (Dunavant): (r, s, w), Σw = 1
const _TRI7 = let a1 = (6 - sqrt(15)) / 21, b1 = (9 + 2sqrt(15)) / 21,
                  a2 = (6 + sqrt(15)) / 21, b2 = (9 - 2sqrt(15)) / 21,
                  w1 = (155 - sqrt(15)) / 1200, w2 = (155 + sqrt(15)) / 1200
    ((1 / 3, 1 / 3, 9 / 40),
     (a1, a1, w1), (b1, a1, w1), (a1, b1, w1),
     (a2, a2, w2), (b2, a2, w2), (a2, b2, w2))
end

# MITC3+ tying points (Lee, Lee & Bathe 2014, Table 1); d = 1/10 000
const _MITC3P_D = 1e-4
const _MITC3P_A = (1 / 6, 2 / 3)
const _MITC3P_B = (2 / 3, 1 / 6)
const _MITC3P_C = (1 / 6, 1 / 6)
const _MITC3P_DP = (1 / 3 + _MITC3P_D, 1 / 3 - 2_MITC3P_D)
const _MITC3P_EP = (1 / 3 - 2_MITC3P_D, 1 / 3 + _MITC3P_D)
const _MITC3P_FP = (1 / 3 + _MITC3P_D, 1 / 3 + _MITC3P_D)

# rotation interpolation f₁…f₄ with the cubic bubble f₄ = 27rs(1−r−s), and
# its natural derivatives (eq 7)
@inline function _mitc3p_f(r, s)
    f4 = 27r * s * (1 - r - s)
    f4r = 27s * (1 - 2r - s)
    f4s = 27r * (1 - r - 2s)
    f = SVector(1 - r - s - f4 / 3, r - f4 / 3, s - f4 / 3, f4)
    fr = SVector(-1 - f4r / 3, 1 - f4r / 3, -f4r / 3, f4r)
    fs = SVector(-1 - f4s / 3, -f4s / 3, 1 - f4s / 3, f4s)
    return f, fr, fs
end

# covariant transverse shear rows e_r = w,r + β·x,r and e_s = w,s + β·x,s
# (engineering) at (r, s) over the 11 plate DOFs (w, θx, θy)×3, βx₄, βy₄;
# x,r = (x21, y21), x,s = (x31, y31) are constant on the flat triangle
@inline function _mitc3p_covariant(r, s, xr, xs)
    f, _, _ = _mitc3p_f(r, s)
    hr = (-1.0, 1.0, 0.0)
    hs = (-1.0, 0.0, 1.0)
    # β·x,d = βx·x_d + βy·y_d with βx = Σ f θy + f₄βx₄ and βy = −Σ f θx + f₄βy₄
    row(hd, t) = vcat(
        SVector(hd[1], -f[1] * t[2], f[1] * t[1]), SVector(hd[2], -f[2] * t[2], f[2] * t[1]),
        SVector(hd[3], -f[3] * t[2], f[3] * t[1]), SVector(f[4] * t[1], f[4] * t[2]))
    return row(hr, xr), row(hs, xs)
end

"""
    mitc3p_plate_stiffness(xy, D, S) -> SMatrix{9,9}

Local MITC3+ plate stiffness of a flat triangle over the DOFs `(w, θx, θy)` of
nodes 1…3 (Lee, Lee & Bathe 2014): Reissner–Mindlin bending with the
rotations enriched by a cubic bubble (two internal DOFs, statically condensed)
and the MITC3+ assumed transverse shear strain field, which removes shear
locking and keeps the element spatially isotropic.

- `xy`: local in-plane node coordinates, counter-clockwise [length]
- `D`: bending rigidity [force·length]
- `S`: transverse shear rigidity [force/length]
"""
mitc3p_plate_stiffness(xy, D::AbstractMatrix, S::AbstractMatrix) =
    _condense(_mitc3p_K11(xy, D, S), SOneTo(9), SVector(10, 11))

# the 11×11 plate stiffness before condensing the bubble (also the object of
# the published single-element eigenvalue test)
function _mitc3p_K11(xy, D::AbstractMatrix, S::AbstractMatrix)
    prep = _mitc3p_prep(xy)
    Dm, Sm = SMatrix{3,3}(D), SMatrix{2,2}(S)
    Ar = prep[2]
    K = zero(SMatrix{11,11,promote_type(eltype(Dm), eltype(Sm), typeof(Ar))})
    for (r, s, w) in _TRI7
        Bb, Bs = _mitc3p_B(prep, r, s)
        K += (w * Ar) * (Bb' * Dm * Bb + Bs' * Sm * Bs)
    end
    return K
end

# per-element geometry of MITC3+: inverse Jacobian, area and the assumed-shear
# field's constant rows and slope row (eq 17)
function _mitc3p_prep(xy)
    (x1, y1), (x2, y2), (x3, y3) = xy[1], xy[2], xy[3]
    xr, xs = SVector(x2 - x1, y2 - y1), SVector(x3 - x1, y3 - y1)
    detJ = xr[1] * xs[2] - xr[2] * xs[1]              # = 2·area
    # Cartesian from natural derivatives: [∂x; ∂y] = J⁻¹[∂r; ∂s], J = [x,r y,r; x,s y,s]
    Ji = SMatrix{2,2}(xs[2], -xs[1], -xr[2], xr[1]) / detJ   # column-major: [y,s −y,r; −x,s x,r]/detJ
    erA, esA = _mitc3p_covariant(_MITC3P_A..., xr, xs)
    erB, esB = _mitc3p_covariant(_MITC3P_B..., xr, xs)
    erC, esC = _mitc3p_covariant(_MITC3P_C..., xr, xs)
    erD, _ = _mitc3p_covariant(_MITC3P_DP..., xr, xs)
    _, esE = _mitc3p_covariant(_MITC3P_EP..., xr, xs)
    erF, esF = _mitc3p_covariant(_MITC3P_FP..., xr, xs)
    ĉ = (erF - erD) - (esF - esE)
    er0 = (2 / 3) * (erB - esB / 2) + (erC + esC) / 3
    es0 = (2 / 3) * (esA - erA / 2) + (erC + esC) / 3
    return Ji, detJ / 2, er0, es0, ĉ
end

# bending (3×11) and assumed transverse shear (2×11) operators at (r, s)
@inline function _mitc3p_B(prep, r, s)
    Ji, _, er0, es0, ĉ = prep
    _, fr, fs = _mitc3p_f(r, s)
    fx = Ji[1, 1] * fr + Ji[1, 2] * fs
    fy = Ji[2, 1] * fr + Ji[2, 2] * fs
    # κ = (βx,x, βy,y, βx,y + βy,x)
    blk(a) = SMatrix{3,3}(0, 0, 0, 0, -fy[a], -fx[a], fx[a], 0, fy[a])
    bub = SMatrix{3,2}(fx[4], 0, fy[4], 0, fy[4], fx[4])
    Bb = hcat(blk(1), blk(2), blk(3), bub)
    er = er0 + ĉ * ((3s - 1) / 3)
    es = es0 + ĉ * ((1 - 3r) / 3)
    Bs = Ji * vcat(er', es')
    return Bb, Bs
end

# ─────────────────────────────────────────────────────────────────────────────
# Geometry and the global 18×18 triangle stiffness
# ─────────────────────────────────────────────────────────────────────────────

"""
    shell_frame(X1, X2, X3) -> (R, xy, warp)

Local frame of a 3-node flat shell (always flat): rows of `R` are `e₁` along
edge 1 (node 1 → node 2), the normal `e₃ = (X2−X1)×(X3−X1)/|…|` and
`e₂ = e₃ × e₁`; `xy` are the in-plane node coordinates relative to the
centroid [length]; `warp` is 0. Counter-clockwise nodes about `e₃` by
construction.
"""
function shell_frame(X1::AbstractVector, X2::AbstractVector, X3::AbstractVector)
    x1, x2, x3 = _sv3(X1), _sv3(X2), _sv3(X3)
    c = (x1 + x2 + x3) / 3
    nrm = cross(x2 - x1, x3 - x1)
    e3 = nrm / norm(nrm)
    v = x2 - x1
    e1 = v / norm(v)
    e2 = cross(e3, e1)
    R = vcat(e1', e2', e3')
    xy = map(d -> SVector(dot(e1, d), dot(e2, d)), (x1 - c, x2 - c, x3 - c))
    return R, xy, zero(eltype(R))
end

# node-major local layout (u,v,w,θx,θy,ψ)×3 ↔ membrane (u,v,ψ)×3 and plate
# (w,θx,θy)×3, padded with a zero 10th row/column (as `_IDX_MEMBRANE` for quads)
const _IDX_MEMBRANE3 = SVector{18,Int}([(1, 2, 10, 10, 10, 3)[(j-1)%6+1] == 10 ? 10 :
    3 * ((j - 1) ÷ 6) + (1, 2, 10, 10, 10, 3)[(j-1)%6+1] for j in 1:18])
const _IDX_PLATE3 = SVector{18,Int}([(10, 10, 1, 2, 3, 10)[(j-1)%6+1] == 10 ? 10 :
    3 * ((j - 1) ÷ 6) + (10, 10, 1, 2, 3, 10)[(j-1)%6+1] for j in 1:18])

@inline function _pad10(K::SMatrix{9,9,T}) where {T}
    z = zero(SVector{9,T})
    return vcat(hcat(K, z), hcat(z', zero(T)))
end

@inline _shell_local3(Km, Kp) = _pad10(Km)[_IDX_MEMBRANE3, _IDX_MEMBRANE3] + _pad10(Kp)[_IDX_PLATE3, _IDX_PLATE3]

# Tᵀ·K·T with T = diag(R, …, R) (6 blocks), blockwise as `_rotate_shell`
@inline function _rotate_shell3(K::SMatrix{18,18}, R::SMatrix{3,3})
    B = reshape(R' * reshape(K, Size(3, 108)), Size(18, 18))
    C = reshape(R' * reshape(B', Size(3, 108)), Size(18, 18))
    return C'
end

# rotate a node-major 18-vector blockwise (local = R·global per 3-block)
@inline _rotate18(R::SMatrix{3,3}, v::SVector{18}) = vec(R * reshape(v, Size(3, 6)))

"""
    shell_stiffness(section, X1, X2, X3) -> SMatrix{18,18}

Global stiffness of a 3-node flat shell (ANDES-OPT membrane ⊕ MITC3+ plate)
given its node positions; local matrices are built in the element plane
([`shell_frame`](@ref)) and rotated blockwise. The drilling rotation is a
genuine membrane DOF of the ANDES triangle, so the section's drilling penalty
`drilling_rigidity` is not used (it belongs to the ITW quad).
"""
shell_stiffness(section::AbstractShellSection, X1, X2, X3) =
    _shell_K_plain(_shell_rigidities(section)..., X1, X2, X3)

# the triangle kernel on plain data (rigidities + positions): the AD rule's target
function _shell_K_plain(A, D, S, γt, X1, X2, X3)
    R, xy, _ = shell_frame(X1, X2, X3)
    Km = andes_membrane_stiffness(xy, A)
    Kp = mitc3p_plate_stiffness(xy, D, S)
    return _rotate_shell3(_shell_local3(Km, Kp), R)
end

_tri_area(xy) = ((xy[2] - xy[1])[1] * (xy[3] - xy[1])[2] - (xy[2] - xy[1])[2] * (xy[3] - xy[1])[1]) / 2
