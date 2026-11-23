#=
Flat-shell element kernels — pure functions of node positions and section
rigidities (shared by the in-place and differentiable assembly paths).

Formulations and every sign convention: docs/shells/DERIVATION.md.
  membrane: Ibrahimbegović–Taylor–Wilson (1990) drilling quad, D-type
            (3×3 Gauss on K with a condensed bubble; 1-point drilling penalty)
  plate:    MITC4 (Dvorkin & Bathe 1984; Bathe & Dvorkin 1985)

Node order: counter-clockwise about the element normal; natural coordinates
(ξ, η) of nodes 1…4 are (−1,−1), (1,−1), (1,1), (−1,1). Edge k joins node k
to node k+1.

Local per-node DOF order (matches the global 6-slot layout so the global
transformation is a plain block rotation): (u, v, w, θx, θy, ψ), where ψ = θz
is the drilling rotation. The membrane uses (u, v, ψ), the plate (w, θx, θy).
=#

const _Q4_ξ = (-1.0, 1.0, 1.0, -1.0)
const _Q4_η = (-1.0, -1.0, 1.0, 1.0)
const _GAUSS3_X = (-sqrt(3 / 5), 0.0, sqrt(3 / 5))
const _GAUSS3_W = (5 / 9, 8 / 9, 5 / 9)
const _GAUSS2_X = (-1 / sqrt(3), 1 / sqrt(3))

# ─────────────────────────────────────────────────────────────────────────────
# Geometry: local frame of a (possibly warped) quadrilateral
# ─────────────────────────────────────────────────────────────────────────────

"""
    shell_frame(X1, X2, X3, X4) -> (R, xy, warp)

Local frame of a 4-node flat shell from its global node positions:

- `R::SMatrix{3,3}`: rows are the local axes `e₁, e₂, e₃` in global
  coordinates (`v_local = R·v_global`). The normal `e₃` is the normalized
  cross product of the diagonals (well defined for warped quads); `e₁` is
  edge 1 (node 1 → node 2) projected onto the mean plane; `e₂ = e₃ × e₁`.
- `xy::NTuple{4,SVector{2}}`: in-plane local coordinates of the nodes
  relative to the centroid (warped nodes projected onto the mean plane)
  [length].
- `warp::T`: largest out-of-plane offset of a node from the mean plane,
  divided by `√(area)` (dimensionless warping measure; 0 for a flat quad).
"""
# identity for static 3-vectors (the AD path): `SVector{3}(::SVector)` has no
# working Zygote pullback, so never convert what is already static
@inline _sv3(x::SVector{3}) = x
@inline _sv3(x::AbstractVector) = SVector{3}(x)

function shell_frame(X1::AbstractVector, X2::AbstractVector, X3::AbstractVector, X4::AbstractVector)
    x1, x2, x3, x4 = _sv3(X1), _sv3(X2), _sv3(X3), _sv3(X4)
    c = (x1 + x2 + x3 + x4) / 4
    nrm = cross(x3 - x1, x4 - x2)
    e3 = nrm / norm(nrm)
    v = x2 - x1
    w = v - dot(v, e3) * e3
    e1 = w / norm(w)
    e2 = cross(e3, e1)
    R = vcat(e1', e2', e3')
    p = (x1 - c, x2 - c, x3 - c, x4 - c)
    xy = map(d -> SVector(dot(e1, d), dot(e2, d)), p)
    area = norm(nrm) / 2
    warp = maximum(d -> abs(dot(e3, d)), p) / sqrt(area)
    return R, xy, warp
end

# ─────────────────────────────────────────────────────────────────────────────
# Shape functions (natural derivatives) at (ξ, η)
# ─────────────────────────────────────────────────────────────────────────────

# corner (bilinear) functions: N, ∂N/∂ξ, ∂N/∂η
@inline function _q4_N(ξ, η)
    N = SVector(ntuple(a -> (1 + _Q4_ξ[a] * ξ) * (1 + _Q4_η[a] * η) / 4, 4))
    Nξ = SVector(ntuple(a -> _Q4_ξ[a] * (1 + _Q4_η[a] * η) / 4, 4))
    Nη = SVector(ntuple(a -> _Q4_η[a] * (1 + _Q4_ξ[a] * ξ) / 4, 4))
    return N, Nξ, Nη
end

# edge (serendipity midside) functions NSₖ, edge k = node k → k+1
@inline function _q4_NS(ξ, η)
    Sξ = SVector(-ξ * (1 - η), (1 - η^2) / 2, -ξ * (1 + η), -(1 - η^2) / 2)
    Sη = SVector(-(1 - ξ^2) / 2, -(1 + ξ) * η, (1 - ξ^2) / 2, -(1 - ξ) * η)
    return Sξ, Sη
end

# Jacobian J = [x,ξ y,ξ; x,η y,η] and the Cartesian-derivative map
@inline function _q4_jacobian(xy, Nξ, Nη)
    xξ = sum(a -> Nξ[a] * xy[a][1], 1:4)
    yξ = sum(a -> Nξ[a] * xy[a][2], 1:4)
    xη = sum(a -> Nη[a] * xy[a][1], 1:4)
    yη = sum(a -> Nη[a] * xy[a][2], 1:4)
    detJ = xξ * yη - yξ * xη
    return xξ, yξ, xη, yη, detJ
end

# (∂f/∂x, ∂f/∂y) from (∂f/∂ξ, ∂f/∂η) for vectors of functions
@inline function _cartesian(fξ, fη, xξ, yξ, xη, yη, detJ)
    fx = (yη .* fξ .- yξ .* fη) ./ detJ
    fy = (xξ .* fη .- xη .* fξ) ./ detJ
    return fx, fy
end

# ─────────────────────────────────────────────────────────────────────────────
# Membrane: ITW drilling quad (D-type)
# ─────────────────────────────────────────────────────────────────────────────

# Allman edge vectors cₖ = lₖ·nₖ = (Δy, −Δx) of edge k (outward for CCW nodes)
@inline _edge_vector(p, q) = SVector(q[2] - p[2], p[1] - q[1])
@inline _q4_edge_vectors(xy) = (_edge_vector(xy[1], xy[2]), _edge_vector(xy[2], xy[3]),
    _edge_vector(xy[3], xy[4]), _edge_vector(xy[4], xy[1]))

# drilling shape-vector derivatives of node a: dₐ = (c_{a−1}·NS_{a−1} − cₐ·NSₐ)/8
@inline function _drill_derivs(c, Sx, Sy, a::Int)
    p = mod1(a - 1, 4)
    dx = (c[p] * Sx[p] - c[a] * Sx[a]) / 8     # ∂dₐ/∂x  (2-vector: u, v parts)
    dy = (c[p] * Sy[p] - c[a] * Sy[a]) / 8     # ∂dₐ/∂y
    return dx, dy
end

# 3×16 strain operator at a point: (u,v,ψ) × 4 nodes, the 2 bubble DOFs, then
# the 2 complementary incompatible modes (Taylor-corrected derivatives Ex, Ey)
@inline function _itw_B(Nx, Ny, Bx, By, dd, Ix, Iy)
    blk(a) = SMatrix{3,3}(Nx[a], zero(Nx[a]), Ny[a],          # column u
        zero(Nx[a]), Ny[a], Nx[a],                             # column v
        dd[a][1][1], dd[a][2][2], dd[a][2][1] + dd[a][1][2])   # column ψ
    bub = SMatrix{3,2}(Bx, zero(Bx), By, zero(Bx), By, Bx)
    inc = SMatrix{3,2}(Ix[1], zero(Bx), Iy[1], zero(Bx), Iy[2], Ix[2])
    return hcat(blk(1), blk(2), blk(3), blk(4), bub, inc)
end

"""
    itw_membrane_stiffness(xy, A, γt; enhanced = true) -> SMatrix{12,12}

Local membrane stiffness of the Ibrahimbegović–Taylor–Wilson drilling quad
(D-type, ITW eq. 3.24) over the DOFs `(u, v, ψ)` of nodes 1…4.

- `xy`: local in-plane node coordinates (counter-clockwise) [length]
- `A`: membrane rigidity [force/length]
- `γt`: drilling penalty rigidity [force/length]

`K = K̂ + P`: `K̂` integrates `[B G]ᵀA[B G]` by 3×3 Gauss with internal modes
statically condensed; `P` integrates the penalty `γt·∫(skew∇u − ψ)²` with a
single point.

Internal modes: ITW's hierarchical bubble `(1−ξ²)(1−η²)`, plus — when
`enhanced` (the default) — two incompatible modes `u ∝ 1−ξ²`, `v ∝ 1−η²`
with Taylor–Beresford–Wilson (1976) derivatives (centre Jacobian, scaled by
`detJ₀/detJ`, so the patch test still passes on distorted meshes). Allman's
edge modes already span the NORMAL edge quadratics; these two add the
complementary tangential ones, which the Poisson strain of in-plane bending
needs: the enhanced element is exact in pure bending on rectangles (the
published element is 2.6 % stiff at ν = 0.3) and more accurate on coarse
wall meshes. `enhanced = false` reproduces the element as published. See
`docs/shells/DERIVATION.md` §2.1 for the benchmark comparison.
"""
function itw_membrane_stiffness(xy, A::AbstractMatrix, γt::Real; enhanced::Bool = true)
    c = _q4_edge_vectors(xy)
    _, N0ξ, N0η = _q4_N(0.0, 0.0)
    j0 = _q4_jacobian(xy, N0ξ, N0η)
    K = _itw_gauss_term(xy, c, j0, A, 1, 1)
    for i in 1:3, j in 1:3
        (i, j) == (1, 1) && continue
        K += _itw_gauss_term(xy, c, j0, A, i, j)
    end
    Kc = enhanced ? _condense(K, SOneTo(12), SVector(13, 14, 15, 16)) :
         _condense(K, SOneTo(12), SVector(13, 14))
    h, detJ0 = _itw_penalty_row(xy, c)
    return Kc + (γt * 4 * detJ0) * (h * h')
end

# static condensation K̂ = Kaa − Kab Kbb⁻¹ Kba of the internal modes `b`
@inline function _condense(K, a, b)
    Kab = K[a, b]
    return K[a, a] - Kab * (K[b, b] \ Kab')
end

# one 3×3-Gauss contribution w·detJ·BᵀAB (16×16, incl. internal modes) — a
# function barrier so the per-point closures capture only local, never-
# reassigned variables (reassigned captures get boxed → allocations)
@inline function _itw_gauss_term(xy, c, j0, A, i::Int, j::Int)
    ξ, η = _GAUSS3_X[i], _GAUSS3_X[j]
    _, Nξ, Nη = _q4_N(ξ, η)
    Sξ, Sη = _q4_NS(ξ, η)
    xξ, yξ, xη, yη, detJ = _q4_jacobian(xy, Nξ, Nη)
    Nx, Ny = _cartesian(Nξ, Nη, xξ, yξ, xη, yη, detJ)
    Sx, Sy = _cartesian(Sξ, Sη, xξ, yξ, xη, yη, detJ)
    Bx, By = _cartesian(-2ξ * (1 - η^2), -2η * (1 - ξ^2), xξ, yξ, xη, yη, detJ)
    # incompatible modes (1−ξ², 1−η²): centre-Jacobian derivatives × detJ₀/detJ
    Ix, Iy = _cartesian(SVector(-2ξ, zero(ξ)), SVector(zero(η), -2η), j0...)
    s = j0[5] / detJ
    dd = (_drill_derivs(c, Sx, Sy, 1), _drill_derivs(c, Sx, Sy, 2),
          _drill_derivs(c, Sx, Sy, 3), _drill_derivs(c, Sx, Sy, 4))
    B = _itw_B(Nx, Ny, Bx, By, dd, s * Ix, s * Iy)
    return (_GAUSS3_W[i] * _GAUSS3_W[j] * detJ) * (B' * A * B)
end

# penalty row h (12) of skew∇u − ψ at the element centre (∇NB = 0 there, and
# the incompatible modes carry no rotation into the penalty), and detJ
@inline function _itw_penalty_row(xy, c)
    N, Nξ, Nη = _q4_N(0.0, 0.0)
    Sξ, Sη = _q4_NS(0.0, 0.0)
    xξ, yξ, xη, yη, detJ = _q4_jacobian(xy, Nξ, Nη)
    Nx, Ny = _cartesian(Nξ, Nη, xξ, yξ, xη, yη, detJ)
    Sx, Sy = _cartesian(Sξ, Sη, xξ, yξ, xη, yη, detJ)
    h = vcat(_itw_h(N, Nx, Ny, c, Sx, Sy, 1), _itw_h(N, Nx, Ny, c, Sx, Sy, 2),
             _itw_h(N, Nx, Ny, c, Sx, Sy, 3), _itw_h(N, Nx, Ny, c, Sx, Sy, 4))
    return h, detJ
end

@inline function _itw_h(N, Nx, Ny, c, Sx, Sy, a::Int)
    dx, dy = _drill_derivs(c, Sx, Sy, a)
    return SVector(-Ny[a] / 2, Nx[a] / 2, (dx[2] - dy[1]) / 2 - N[a])
end

# ─────────────────────────────────────────────────────────────────────────────
# Assembly of the local 24×24 and rotation to global
# ─────────────────────────────────────────────────────────────────────────────

# index maps of the node-major local layout (u,v,w,θx,θy,ψ)×4 into the
# membrane (u,v,ψ)×4 and plate (w,θx,θy)×4 matrices, padded with a zero 13th
# row/column that absorbs the other part's slots — embedding is then one
# static gather per part instead of a 0/1 matrix product
const _IDX_MEMBRANE = SVector{24,Int}([(1, 2, 13, 13, 13, 3)[(j-1)%6+1] == 13 ? 13 :
    3 * ((j - 1) ÷ 6) + (1, 2, 13, 13, 13, 3)[(j-1)%6+1] for j in 1:24])
const _IDX_PLATE = SVector{24,Int}([(13, 13, 1, 2, 3, 13)[(j-1)%6+1] == 13 ? 13 :
    3 * ((j - 1) ÷ 6) + (13, 13, 1, 2, 3, 13)[(j-1)%6+1] for j in 1:24])

@inline function _pad13(K::SMatrix{12,12,T}) where {T}
    z = zero(SVector{12,T})
    return vcat(hcat(K, z), hcat(z', zero(T)))
end

# local 24×24 from the membrane and plate 12×12s
@inline _shell_local(Km, Kp) = _pad13(Km)[_IDX_MEMBRANE, _IDX_MEMBRANE] + _pad13(Kp)[_IDX_PLATE, _IDX_PLATE]

# Tᵀ·K·T with T = diag(R, …, R) (8 blocks), done blockwise: a column-major
# reshape of the 24×24 into 3×192 stacks every 3-row block as a column, so
# one 3×3 product rotates all row blocks at once (Tᵀ·K); K symmetric ⇒ the
# column rotation is the same operation on the transpose
@inline function _rotate_shell(K::SMatrix{24,24}, R::SMatrix{3,3})
    B = reshape(R' * reshape(K, Size(3, 192)), Size(24, 24))      # Tᵀ K
    C = reshape(R' * reshape(B', Size(3, 192)), Size(24, 24))     # Tᵀ (Tᵀ K)ᵀ = (Tᵀ K T)ᵀ
    return C'
end

"""
    shell_stiffness(section, X1, X2, X3, X4) -> SMatrix{24,24}

Global stiffness of a 4-node flat shell (membrane ⊕ plate) given its node
positions. Local matrices are built in the element plane
([`shell_frame`](@ref)) and rotated blockwise: `K = Tᵀ K_local T`,
`T = diag(R, …, R)`.
"""
function shell_stiffness(section::AbstractShellSection, X1, X2, X3, X4)
    R, xy, _ = shell_frame(X1, X2, X3, X4)
    Km = itw_membrane_stiffness(xy, membrane_rigidity(section), drilling_rigidity(section))
    Kp = _plate_stiffness(xy, section)
    return _rotate_shell(_shell_local(Km, Kp), R)
end

# plate (w, θx, θy)×4 stiffness — MITC4 lands in S3; until then the element is
# membrane-only and its out-of-plane DOFs carry no stiffness
_plate_stiffness(xy, section) = zero(SMatrix{12,12,promote_type(eltype(xy[1]), eltype(membrane_rigidity(section)))})
