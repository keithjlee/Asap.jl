#=
Shell element loads (consistent nodal loads) and stress-resultant recovery —
pure kernels of the local geometry, rigidities and local displacements.
See docs/shells/DERIVATION.md §2 (membrane interpolation) and §4.
=#

# the 16×16 membrane matrix BEFORE condensation: (u,v,ψ)×4, bubble (2),
# incompatible modes (2)
function _itw_K16(xy, A::AbstractMatrix)
    c = _q4_edge_vectors(xy)
    _, N0ξ, N0η = _q4_N(0.0, 0.0)
    j0 = _q4_jacobian(xy, N0ξ, N0η)
    K = _itw_gauss_term(xy, c, j0, A, 1, 1)
    for i in 1:3, j in 1:3
        (i, j) == (1, 1) && continue
        K += _itw_gauss_term(xy, c, j0, A, i, j)
    end
    return K
end

# 3×16 membrane strain operator at (ξ, η) — the same B the stiffness integrates
function _itw_B16(xy, ξ, η)
    c = _q4_edge_vectors(xy)
    _, N0ξ, N0η = _q4_N(0.0, 0.0)
    j0 = _q4_jacobian(xy, N0ξ, N0η)
    _, Nξ, Nη = _q4_N(ξ, η)
    Sξ, Sη = _q4_NS(ξ, η)
    xξ, yξ, xη, yη, detJ = _q4_jacobian(xy, Nξ, Nη)
    Nx, Ny = _cartesian(Nξ, Nη, xξ, yξ, xη, yη, detJ)
    Sx, Sy = _cartesian(Sξ, Sη, xξ, yξ, xη, yη, detJ)
    Bx, By = _cartesian(-2ξ * (1 - η^2), -2η * (1 - ξ^2), xξ, yξ, xη, yη, detJ)
    Ix, Iy = _cartesian(SVector(-2ξ, zero(ξ)), SVector(zero(η), -2η), j0...)
    s = j0[5] / detJ
    dd = (_drill_derivs(c, Sx, Sy, 1), _drill_derivs(c, Sx, Sy, 2),
          _drill_derivs(c, Sx, Sy, 3), _drill_derivs(c, Sx, Sy, 4))
    return _itw_B(Nx, Ny, Bx, By, dd, s * Ix, s * Iy)
end

# serendipity midside function VALUES NSₖ(ξ, η), edge k = node k → k+1
@inline _q4_NSval(ξ, η) = SVector((1 - ξ^2) * (1 - η) / 2, (1 + ξ) * (1 - η^2) / 2,
    (1 - ξ^2) * (1 + η) / 2, (1 - ξ) * (1 - η^2) / 2)

"""
    shell_consistent_load(xy, q) -> (f24, fb)

Consistent local nodal loads of a uniform surface traction `q = (qx, qy, qz)`
[force/length²] in the element's local axes, from the element's own
interpolation (3×3 Gauss):

- membrane: `∫Nₐ q∥` on (u, v), and `∫dₐ·q∥` on the drilling DOF ψ (the
  Allman edge modes do work under in-plane load);
- plate: `∫Nₐ qz` on w (MITC4 rotations are independent fields — a uniform
  pressure does no work on them);
- `fb`: the load on the internal modes (bubble `∫NB q∥`; the incompatible
  modes are strain enhancements and carry none), used for condensation.

Returns `f24` in the node-major local layout `(u, v, w, θx, θy, ψ)×4` BEFORE
condensation of the internal modes, and `fb` (4: bubble x, y, incompatible 1, 2).
"""
function shell_consistent_load(xy, q::SVector{3})
    c = _q4_edge_vectors(xy)
    T = promote_type(eltype(xy[1]), eltype(q))
    fu = zero(SVector{4,T}); fv = zero(SVector{4,T}); fψ = zero(SVector{4,T}); fw = zero(SVector{4,T})
    fbx = zero(T); fby = zero(T)
    for i in 1:3, j in 1:3
        ξ, η = _GAUSS3_X[i], _GAUSS3_X[j]
        N, Nξ, Nη = _q4_N(ξ, η)
        detJ = _q4_jacobian(xy, Nξ, Nη)[5]
        S = _q4_NSval(ξ, η)
        w = _GAUSS3_W[i] * _GAUSS3_W[j] * detJ
        fu += (w * q[1]) * N
        fv += (w * q[2]) * N
        fw += (w * q[3]) * N
        # dₐ·q∥ with dₐ = (c_{a−1} NS_{a−1} − cₐ NSₐ)/8
        cq = SVector(dot(c[1], q[SOneTo(2)]), dot(c[2], q[SOneTo(2)]),
                     dot(c[3], q[SOneTo(2)]), dot(c[4], q[SOneTo(2)]))
        fψ += (w / 8) * SVector(cq[4] * S[4] - cq[1] * S[1], cq[1] * S[1] - cq[2] * S[2],
                                cq[2] * S[2] - cq[3] * S[3], cq[3] * S[3] - cq[4] * S[4])
        nb = (1 - ξ^2) * (1 - η^2)
        fbx += w * nb * q[1]
        fby += w * nb * q[2]
    end
    z = zero(T)
    f24 = SVector{24,T}(fu[1], fv[1], fw[1], z, z, fψ[1], fu[2], fv[2], fw[2], z, z, fψ[2],
                        fu[3], fv[3], fw[3], z, z, fψ[3], fu[4], fv[4], fw[4], z, z, fψ[4])
    return f24, SVector{4,T}(fbx, fby, z, z)
end

const _IDX_M12 = SVector(1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20, 24)   # (u,v,ψ)×4 in the 24-layout

"""
    shell_load_vector(section, X1, X2, X3, X4, q_global) -> (q̃_local, R)

The shell's fixed-end-force vector (Asap convention: `Pf` and the `q` in the
force recovery `f = K·u + q̃`) for a uniform surface traction given in GLOBAL
axes [force/length²]: the NEGATED consistent load, with the internal-mode
load condensed onto the nodes, `f̂ = fₐ − K_ab K_bb⁻¹ f_b`. Local
node-major 24-layout; `R` is the element frame (rows = local axes).
"""
function shell_load_vector(section::AbstractShellSection, X1, X2, X3, X4, q::SVector{3})
    R, xy, _ = shell_frame(X1, X2, X3, X4)
    f24, fb = shell_consistent_load(xy, R * q)
    K16 = _itw_K16(xy, membrane_rigidity(section))
    b = SVector(13, 14, 15, 16)
    corr = K16[_SONE12, b] * (K16[b, b] \ fb)          # K_ab K_bb⁻¹ f_b on (u,v,ψ)×4
    Δ = vcat(corr, zero(eltype(corr)))[_IDX_MEMBRANE]   # scatter (u,v,ψ)×4 into the 24-layout
    return -(f24 - Δ), R
end

const _SONE12 = SOneTo(12)

# rotate a node-major 24-vector blockwise by R (local = R·global per 3-block)
@inline _rotate24(R::SMatrix{3,3}, v::SVector{24}) =
    vec(R * reshape(v, Size(3, 8)))

"""
    shell_resultants_local(section, xy, a24, fb; at = (0.0, 0.0)) -> NamedTuple

Stress resultants at natural point `at` from LOCAL nodal displacements `a24`
(node-major `(u, v, w, θx, θy, ψ)×4`) and the internal-mode load `fb`:

- `N = (Nxx, Nyy, Nxy)` [force/length] — membrane forces, tension positive
- `M = (Mxx, Myy, Mxy)` [force] — moments per unit length; positive `Mxx`
  puts the `+e₃` face in tension along `e₁`
- `Q = (Qx, Qy)` [force/length] — transverse shears (MITC4 field)

The membrane's condensed modes are recovered as
`α = K_bb⁻¹ (f_b − K_ba a)` so `N` is the field the stiffness integrated.
"""
function shell_resultants_local(section::AbstractShellSection, xy, a24::SVector{24}, fb::SVector{4};
    at = (0.0, 0.0))
    ξ, η = at
    A = membrane_rigidity(section)
    am = a24[_IDX_M12]
    ap = a24[SVector(3, 4, 5, 9, 10, 11, 15, 16, 17, 21, 22, 23)]
    K16 = _itw_K16(xy, A)
    b = SVector(13, 14, 15, 16)
    α = K16[b, b] \ (fb - K16[b, _SONE12] * am)
    ε = _itw_B16(xy, ξ, η) * vcat(am, α)
    N = A * ε
    # plate: curvature and MITC shear at (ξ, η)
    _, Nξ, Nη = _q4_N(ξ, η)
    xξ, yξ, xη, yη, detJ = _q4_jacobian(xy, Nξ, Nη)
    Nx, Ny = _cartesian(Nξ, Nη, xξ, yξ, xη, yη, detJ)
    M = bending_rigidity(section) * (_mitc_Bb(Nx, Ny) * ap)
    Q = shear_rigidity(section) * _mitc_gamma(xy, ap, ξ, η)
    return (N = N, M = M, Q = Q)
end

"""
    _mitc_gamma(xy, ap, ξ, η) -> SVector{2}

The MITC4 assumed transverse shear strains `(γxz, γyz)` at `(ξ, η)` from the
plate DOFs `ap = (w, θx, θy)×4` — the same field the stiffness integrates
(covariant strains tied at the edge midpoints, mapped by `J⁻¹`). Verified
against Bathe & Dvorkin (1985) eqs. 8–10 in `test_shell_mitc4_bd85.jl`.
"""
function _mitc_gamma(xy, ap, ξ, η)
    _, Nξ, Nη = _q4_N(ξ, η)
    xξ, yξ, xη, yη, detJ = _q4_jacobian(xy, Nξ, Nη)
    eA = _mitc_tying_row(xy, 0.0, -1.0, 1); eC = _mitc_tying_row(xy, 0.0, 1.0, 1)
    eB = _mitc_tying_row(xy, 1.0, 0.0, 2); eD = _mitc_tying_row(xy, -1.0, 0.0, 2)
    eξ = dot(((1 - η) / 2) * eA + ((1 + η) / 2) * eC, ap)
    eη = dot(((1 + ξ) / 2) * eB + ((1 - ξ) / 2) * eD, ap)
    return SVector((yη * eξ - yξ * eη) / detJ, (xξ * eη - xη * eξ) / detJ)
end
