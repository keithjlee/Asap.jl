#=
Fixed-end force kernels.

One quadrature engine replaces the closed-form fixed-end-force catalog:
the consistent nodal load for ANY distributed intensity is the integral of
the load against the element's shape functions, and for piecewise-linear
intensities that integrand is polynomial of degree ≤ 4 — so 3-point
Gauss-Legendre (exact to degree 5) evaluates it EXACTLY, not approximately.

The engine produces the clamped-clamped local 12-vector; end-condition
condensation (releases / semi-rigid springs) is applied afterwards,
generically, by `condense_fef` — loads know nothing about end conditions.

Sign convention (legacy-compatible): fixed-end forces are the reactions the
clamped supports must supply, i.e. q = −∫ N(x)ᵀ w(x) dx; downstream,
Pf accumulates Λᵀq and the solve uses F = P − Pf, and local force recovery
is f = k·u_local + q̃.
=#

# 3-point Gauss–Legendre nodes/weights on [−1, 1]: exact for degree ≤ 5
const GAUSS3_X = SVector(-sqrt(3 / 5), 0.0, sqrt(3 / 5))
const GAUSS3_W = SVector(5 / 9, 8 / 9, 5 / 9)

"""
    shape_functions(x, L) -> (Na, Nb, H1, H2, H3, H4)

Element shape functions evaluated at local coordinate `x ∈ [0, L]`:
linear (`Na`, `Nb`) for axial/torsional interpolation and cubic Hermite
(`H1`–`H4`) for transverse displacement — `H1`/`H3` pair with end
translations, `H2`/`H4` with end rotations. These are the same functions the
element stiffness is derived from, which is what makes quadrature against
them the *consistent* load lowering.
"""
@inline function shape_functions(x::Real, L::Real)
    ξ = x / L
    Na = 1 - ξ
    Nb = ξ
    H1 = 1 - 3ξ^2 + 2ξ^3
    H2 = L * (ξ - 2ξ^2 + ξ^3)
    H3 = 3ξ^2 - 2ξ^3
    H4 = L * (ξ^3 - ξ^2)
    return Na, Nb, H1, H2, H3, H4
end

"""
    _fef_contribution(c, d, x, L) -> SVector{12}

The local 12-vector contribution of a force density `c·d` (`d` a local
direction 3-vector) acting at station `x ∈ [0, L]`: axial against the linear
shape functions, transverse against the cubic Hermite ones (local DOF order;
y couples +H2/+H4 rotations, z couples −H2/−H4 — the legacy x–z plane sign
convention).

Returned as a fresh immutable `SVector` rather than accumulated into a
scratch buffer: the FEF kernels are on the differentiable pure path
(design-dependent loads), and an in-place accumulation is opaque to
reverse-mode AD (Zygote). Summation order is unchanged, so the stateful
path's numbers are too.
"""
@inline function _fef_contribution(c::Real, d::AbstractVector, x::Real, L::Real)
    Na, Nb, H1, H2, H3, H4 = shape_functions(x, L)
    return SVector(
        c * d[1] * Na, c * d[2] * H1, c * d[3] * H1, zero(c * d[1] * Na),
        -(c * d[3] * H2), c * d[2] * H2,
        c * d[1] * Nb, c * d[2] * H3, c * d[3] * H3, zero(c * d[1] * Nb),
        -(c * d[3] * H4), c * d[2] * H4,
    )
end

"""
    _distributed_fef(t, w, d, L) -> SVector{12}

Clamped-clamped local FEF of a piecewise-linear distributed load: intensity
`w[s] → w[s+1]` [force/length] over each fraction interval `t[s] → t[s+1]` of
the CURRENT length `L` (so the total load scales with `L` — per-unit-length
semantics), along the local direction `d`. Three-point Gauss–Legendre per
segment, exact for the degree-≤4 integrand.
"""
function _distributed_fef(t::AbstractVector, w::AbstractVector, d::AbstractVector, L::Real)
    TQ = promote_type(eltype(w), eltype(d), typeof(L))
    q = zero(SVector{12,TQ})
    for s in 1:length(t)-1
        # degenerate (zero-length) segment: decided on the fractions, which
        # are constants — never on a possibly-Dual length
        t[s+1] - t[s] < eps(float(eltype(t))) && continue
        a = t[s] * L
        b = t[s+1] * L
        h = (b - a) / 2
        m = (a + b) / 2
        for k in 1:3
            x = m + h * GAUSS3_X[k]
            wi = w[s] + (w[s+1] - w[s]) * (x - a) / (b - a)
            c = -GAUSS3_W[k] * h * wi          # FEF = −consistent nodal load
            q = q + _fef_contribution(c, d, x, L)
        end
    end
    return q
end

# direction of an element load in LOCAL coordinates
@inline _local_direction(direction::SVector{3}, coords::Symbol, Λ::SMatrix{3,3}) =
    coords === :local ? direction : Λ * direction

"""
    fixed_end_forces(load, section, ends, x1, x2, rollangle) -> SVector{12}

Clamped-clamped fixed-end force vector in element LOCAL coordinates for an
element load. THE extension point for new element load types: implement
exactly this method; condensation for releases/semi-rigid ends and rotation
to global coordinates are applied generically by the load assembler.

`section`, `x1`, `x2`, `rollangle` describe the element (positions passed explicitly
so the pure AD path can differentiate through geometry); `ends` is accepted
for interface uniformity but the returned vector is always the CLAMPED one.

# Semantics under design change (the pure path re-evaluates FEFs per state)
- distributed loads are per-unit-length intensities over fractions of the
  member: the total load is recomputed from the current length
- point forces/moments sit at a fixed FRACTION of the member, so they ride
  with it
- `SelfWeight` follows the current `ρA` and length
- `:local` directions rotate with the member; `:global` ones do not

All methods are non-mutating (immutable `SVector` arithmetic), so they are
differentiable by ForwardDiff, Zygote and the other AD backends.
"""
function fixed_end_forces(load::DistributedLoad{T}, section::AbstractSection,
    ends::EndConditions, x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real},
    rollangle::Real) where {T}

    L = element_length(x1, x2)
    Λ = local_frame(x1, x2, rollangle)
    d = _local_direction(load.direction, load.coords, Λ)
    return _distributed_fef(load.t, load.w, d, L)
end

function fixed_end_forces(load::PointLoad{T}, section::AbstractSection,
    ends::EndConditions, x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real},
    rollangle::Real) where {T}

    L = element_length(x1, x2)
    Λ = local_frame(x1, x2, rollangle)
    p = load.coords === :local ? load.value : Λ * load.value
    # a point force is the c = −1 limit of the quadrature: one kernel evaluation
    return _fef_contribution(-one(eltype(p)), p, load.position * L, L)
end

"""
    shape_function_slopes(x, L) -> (Na′, Nb′, H1′, H2′, H3′, H4′)

Derivatives (w.r.t. the local coordinate) of [`shape_functions`](@ref) —
the kernel for concentrated MOMENTS, which do work through the rotation
(slope) field rather than the displacement field.
"""
@inline function shape_function_slopes(x::Real, L::Real)
    ξ = x / L
    Na′ = -1 / L
    Nb′ = 1 / L
    H1′ = (-6ξ + 6ξ^2) / L
    H2′ = 1 - 4ξ + 3ξ^2
    H3′ = (6ξ - 6ξ^2) / L
    H4′ = 3ξ^2 - 2ξ
    return Na′, Nb′, H1′, H2′, H3′, H4′
end

function fixed_end_forces(load::PointMoment{T}, section::AbstractSection,
    ends::EndConditions, x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real},
    rollangle::Real) where {T}

    L = element_length(x1, x2)
    Λ = local_frame(x1, x2, rollangle)
    m = load.coords === :local ? load.value : Λ * load.value
    x = load.position * L

    Na, Nb, _ = shape_functions(x, L)
    _, _, H1′, H2′, H3′, H4′ = shape_function_slopes(x, L)

    # slots: torsion (moment about local x) interpolates linearly, like axial;
    # moment about local z pairs with the x–y plane displacement slope; moment
    # about local y is the x–z plane with the flipped rotation pairing
    return SVector(
        zero(m[1] * Na), -m[3] * H1′, m[2] * H1′, -m[1] * Na, -m[2] * H2′, -m[3] * H2′,
        zero(m[1] * Nb), -m[3] * H3′, m[2] * H3′, -m[1] * Nb, -m[2] * H4′, -m[3] * H4′,
    )
end

function fixed_end_forces(load::SelfWeight{T}, section::AbstractSection,
    ends::EndConditions, x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real},
    rollangle::Real) where {T}

    gmag = norm(load.g)
    iszero(gmag) && return zero(SVector{12,promote_type(T, eltype(x1), eltype(x2), typeof(ρA(section)))})
    # a uniform global line load of intensity ρA·|g|·factor along ĝ — lowered
    # through the kernel directly (no DistributedLoad{T} construction, which
    # could not hold a Dual ρA)
    w = ρA(section) * gmag * load.factor
    L = element_length(x1, x2)
    Λ = local_frame(x1, x2, rollangle)
    d = Λ * (load.g ./ gmag)
    return _distributed_fef(SVector(zero(T), one(T)), SVector(w, w), d, L)
end

"""
    condense_fef(q, section, L, ends::EndConditions) -> SVector{12}

Transform a clamped-clamped local fixed-end force vector `q` into the
fixed-end forces consistent with the element's actual end conditions.

Physics: the beam's end rotations connect to the nodes through the end
springs. Holding all nodal DOFs fixed, the beam-side rotations `θb` of each
bending plane relax to

    θb = −(Kθθ + diag(k))⁻¹ qθ

(where `Kθθ` is the clamped bending block's rotation–rotation part and `k`
the end spring pair), and the nodal fixed-end forces become

    q̃θ = qθ + Kθθ·θb          (rotations — zero at an ideal release)
    q̃v = qv + Kvθ·θb          (translations pick up the redistributed shear)

Rigid springs (`Inf`) pass the clamped values through unchanged; zero
springs reproduce the legacy release corrections exactly (characterization-
tested); finite springs give the consistent semi-rigid redistribution.
Axial and torsional FEF components pass through unchanged — no distributed
axial-spring/torque load type exists yet (revisit when one does).
"""
function condense_fef(q::SVector{12,T}, section::AbstractSection, L::Real,
    ends::EndConditions) where {T}

    # x–y plane: (v1, θz1, v2, θz2) = slots (2, 6, 8, 12), standard signs, EIx
    qz = _condense_plane(SVector(q[2], q[6], q[8], q[12]),
        EIx(section), L, ends.e1.kz, ends.e2.kz, false)
    # x–z plane: (v1, θy1, v2, θy2) = slots (3, 5, 9, 11), flipped signs, EIy
    qy = _condense_plane(SVector(q[3], q[5], q[9], q[11]),
        EIy(section), L, ends.e1.ky, ends.e2.ky, true)

    # eltype promotes (not fixed to q's): EI or the springs may carry AD
    # partials that the clamped q (geometry/intensity only) does not
    return SVector(promote(
        q[1], qz[1], qy[1], q[4], qy[2], qz[2],
        q[7], qz[3], qy[3], q[10], qy[4], qz[4],
    )...)
end

function _condense_plane(qp::SVector{4,T}, EI::Real, L::Real,
    k1::Real, k2::Real, flip::Bool) where {T}

    (isinf(k1) && isinf(k2)) && return qp        # clamped: nothing to relax

    Kb = bending_block(EI, L, one(T), one(T); flip=flip)   # clamped block
    Kθθ11, Kθθ12, Kθθ22 = Kb[2, 2], Kb[2, 4], Kb[4, 4]
    qθ1, qθ2 = qp[2], qp[4]

    # solve (Kθθ + diag(k)) θb = −qθ with per-entry Inf handling
    if isinf(k1)
        θ1 = zero(T)
        θ2 = -qθ2 / (Kθθ22 + k2)
    elseif isinf(k2)
        θ1 = -qθ1 / (Kθθ11 + k1)
        θ2 = zero(T)
    else
        a11 = Kθθ11 + k1
        a22 = Kθθ22 + k2
        det = a11 * a22 - Kθθ12^2
        θ1 = (-qθ1 * a22 + qθ2 * Kθθ12) / det
        θ2 = (-qθ2 * a11 + qθ1 * Kθθ12) / det
    end

    return SVector(
        qp[1] + Kb[1, 2] * θ1 + Kb[1, 4] * θ2,
        qθ1 + Kθθ11 * θ1 + Kθθ12 * θ2,
        qp[3] + Kb[3, 2] * θ1 + Kb[3, 4] * θ2,
        qθ2 + Kθθ12 * θ1 + Kθθ22 * θ2,
    )
end
