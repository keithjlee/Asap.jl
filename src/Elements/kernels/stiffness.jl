#=
Element stiffness kernels.

One general formulation — the Monforton–Wu fixity-factor beam with end
springs — replaces the five closed-form release matrices of the legacy
library. The classical releases are exact limits (verified against the
pinned legacy matrices in the characterization suite); intermediate spring
values model semi-rigid connections and are differentiable design variables.

All kernels are pure functions of plain values returning SMatrix — shared by
the in-place assembler and the pure AD path.

Local DOF order (12): [u1x u1y u1z θ1x θ1y θ1z u2x u2y u2z θ2x θ2y θ2z]
Bending in the local x–y plane (EIx, the strong axis) couples DOFs
(u_y, θ_z); bending in the local x–z plane (EIy) couples (u_z, θ_y) with
opposite-signed off-diagonal terms — the legacy sign convention, preserved.
=#

"""
    fixity_factor(k, EI, L) -> p ∈ [0, 1]

Monforton–Wu fixity factor for a rotational end spring of stiffness `k`
[force·length/rad] on a member of flexural rigidity `EI` and length `L`:

    p = 1 / (1 + 3EI / (k·L))

`p = 1` is a rigid connection (`k = ∞`), `p = 0` an ideal hinge (`k = 0`).
The infinite and zero cases branch explicitly so automatic differentiation
never propagates through `Inf` (which would produce NaN partials).
"""
function fixity_factor(k::Real, EI::Real, L::Real)
    isinf(k) && return one(promote_type(typeof(EI), typeof(L)))
    iszero(k) && return zero(promote_type(typeof(EI), typeof(L)))
    return 1 / (1 + 3 * EI / (k * L))
end

"""
    series_rigidity(direct, k1, k2) -> k_eff

Effective stiffness of a member stiffness `direct` (e.g. `EA/L` or `GJ/L`)
in series with two end springs `k1`, `k2` — the axial/torsional analogue of
the fixity factor:

    k_eff = 1 / (1/direct + 1/k1 + 1/k2)

Rigid springs (`Inf`) drop out; a zero spring at either end releases the
action entirely (`k_eff = 0`). Zero/infinite cases branch explicitly for AD
safety.
"""
function series_rigidity(direct::T, k1::Real, k2::Real) where {T<:Real}
    (iszero(k1) || iszero(k2)) && return zero(T)
    c = 1 / direct
    isinf(k1) || (c += 1 / k1)
    isinf(k2) || (c += 1 / k2)
    return 1 / c
end

"""
    bending_block(EI, L, p1, p2; flip = false) -> SMatrix{4,4}

Monforton–Wu condensed bending stiffness on the DOF set
`(v₁, θ₁, v₂, θ₂)` — transverse displacement and bending rotation at each
end — for a beam of flexural rigidity `EI`, length `L`, and end fixity
factors `p1`, `p2` (see [`fixity_factor`](@ref)):

    EI / (L³(4 − p₁p₂)) ×
    ⎡ 12(p₁+p₂+p₁p₂)   6L·p₁(2+p₂)    −12(p₁+p₂+p₁p₂)   6L·p₂(2+p₁) ⎤
    ⎢                  12L²·p₁        −6L·p₁(2+p₂)      6L²·p₁p₂    ⎥
    ⎢       (sym)                     12(p₁+p₂+p₁p₂)    −6L·p₂(2+p₁)⎥
    ⎣                                                   12L²·p₂     ⎦

At `p₁ = p₂ = 1` this is the classical `[12, 6L, 4L², 2L²]·EI/L³` block; at
`p = 0` the corresponding rotation decouples exactly as in the legacy
release matrices.

`flip = true` negates the displacement–rotation coupling terms — the sign
convention for bending in the local x–z plane (positive rotation about
local y pairs with negative transverse-z shear coupling), matching the
legacy `Iy` block.
"""
function bending_block(EI::Real, L::Real, p1::Real, p2::Real; flip::Bool=false)
    den = L^3 * (4 - p1 * p2)
    c = EI / den

    k11 = 12 * (p1 + p2 + p1 * p2) * c
    k12 = 6 * L * p1 * (2 + p2) * c
    k14 = 6 * L * p2 * (2 + p1) * c
    k22 = 12 * L^2 * p1 * c
    k24 = 6 * L^2 * p1 * p2 * c
    k44 = 12 * L^2 * p2 * c

    s = flip ? -1 : 1
    k12 *= s
    k14 *= s
    # (2,3) and (3,4) couplings are −k12, −k14 before the flip; sign carried through

    return SMatrix{4,4}(
        k11, k12, -k11, k14,
        k12, k22, -k12, k24,
        -k11, -k12, k11, -k14,
        k14, k24, -k14, k44,
    )
end

"""
    local_stiffness(section, L, ends::EndConditions) -> SMatrix{12,12}

Element stiffness matrix in LOCAL coordinates for a 3D frame element of
length `L` with cross-section `section` and end conditions `ends`.

Composed from four independent actions, each consuming one rigidity
accessor of the section contract:

- axial (`EA`): series combination with the end axial springs
- torsion (`GJ`): series combination with the end torsional springs
- bending in the local x–y plane (`EIx`): Monforton–Wu block with fixity
  factors from the end `kz` rotational springs
- bending in the local x–z plane (`EIy`): Monforton–Wu block (sign-flipped
  couplings) with fixity factors from the end `ky` springs

The five classical releases reproduce the legacy closed-form matrices
exactly (characterization-tested); any other spring values give a
semi-rigid member.
"""
function local_stiffness(section::AbstractSection, L::Real, ends::EndConditions)
    ea = EA(section)
    eix = EIx(section)
    eiy = EIy(section)
    gj = GJ(section)

    ka = series_rigidity(ea / L, ends.e1.kx, ends.e2.kx)
    kt = series_rigidity(gj / L, ends.e1.kt, ends.e2.kt)

    # x–y plane bending (strong axis EIx): rotations about local z
    pz1 = fixity_factor(ends.e1.kz, eix, L)
    pz2 = fixity_factor(ends.e2.kz, eix, L)
    Bz = bending_block(eix, L, pz1, pz2)

    # x–z plane bending (weak axis EIy): rotations about local y, flipped signs
    py1 = fixity_factor(ends.e1.ky, eiy, L)
    py2 = fixity_factor(ends.e2.ky, eiy, L)
    By = bending_block(eiy, L, py1, py2; flip=true)

    # compose by embedding each action block through constant selector
    # matrices — a pure (mutation-free) expression, so the identical kernel
    # is differentiable by Zygote and fast on the in-place path
    pair = SMatrix{2,2}(1, -1, -1, 1)
    return _SEL_AXIAL * (ka * pair) * _SEL_AXIAL' +
           _SEL_TORSION * (kt * pair) * _SEL_TORSION' +
           _SEL_BEND_Z * Bz * _SEL_BEND_Z' +
           _SEL_BEND_Y * By * _SEL_BEND_Y'
end

# constant selector (embedding) matrices: column j has a 1 at the local DOF
# slot the block's j-th entry maps to
function _selector(idx::NTuple{N,Int}) where {N}
    return SMatrix{12,N,Float64}(ntuple(k -> begin
            i = (k - 1) % 12 + 1
            j = (k - 1) ÷ 12 + 1
            idx[j] == i ? 1.0 : 0.0
        end, 12N))
end

const _SEL_AXIAL = _selector((1, 7))
const _SEL_TORSION = _selector((4, 10))
const _SEL_BEND_Z = _selector((2, 6, 8, 12))    # x–y plane: (v, θz) pairs
const _SEL_BEND_Y = _selector((3, 5, 9, 11))    # x–z plane: (v, θy) pairs
const _SEL_TRANSLATIONS = _selector((1, 2, 3, 7, 8, 9))

"""
    truss_stiffness(section, x1, x2) -> SMatrix{6,6}

GLOBAL-coordinate stiffness of an axial-only (truss) element between
positions `x1` and `x2`: `EA/L · [nnᵀ −nnᵀ; −nnᵀ nnᵀ]` where `n` is the
unit vector along the element.

Returned directly in global coordinates on the 6 translational DOFs
(3 per node) — a truss element never touches rotational DOF slots, which is
what keeps those slots inactive (and out of the solve) in truss-only models.
"""
function truss_stiffness(section::AbstractSection, x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real})
    v = SVector{3}(x2) - SVector{3}(x1)
    L = norm(v)
    n = v / L
    B = (EA(section) / L) * (n * n')
    return [B -B; -B B]
end

"""
    frame_stiffness(section, ends, x1, x2, rollangle) -> SMatrix{12,12}

GLOBAL-coordinate stiffness of a 3D frame element: the local
[`local_stiffness`](@ref) rotated through [`local_frame`](@ref) via
[`transform_to_global`](@ref). Pure function of positions and properties —
the single kernel shared by the in-place and AD assembly paths.
"""
function frame_stiffness(section::AbstractSection, ends::EndConditions,
    x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real}, rollangle::Real)
    L = element_length(x1, x2)
    Λ = local_frame(x1, x2, rollangle)
    return transform_to_global(local_stiffness(section, L, ends), Λ)
end
