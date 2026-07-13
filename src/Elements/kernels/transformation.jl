#=
Geometric transformation kernels.

Pure functions of plain values (positions in, SMatrix out) so the identical
code serves the in-place fast path, the pure AD path, and any scalar type.
=#

"""
    local_frame(x1, x2, Ψ; tol = 1e-6) -> SMatrix{3,3}

Rotation matrix `Λ` from global to element-local coordinates for a line
element from position `x1` to `x2` with roll angle `Ψ`.

The rows of `Λ` are the element's local unit axes expressed in global
coordinates:

- row 1: local x — along the element, from start to end
- row 2: local y — a transverse axis; with `Ψ = π/2` (the constructor
  default) it aligns with "strong-axis bending carries vertical load" for
  typical members
- row 3: local z — completes the right-handed triad

`Ψ` [rad] rotates the local y–z pair about the element axis — the "roll" of
the section. `tol` triggers the special-case frame for members parallel to
the global Y axis, where the general formula degenerates (its denominator
`√(CXx² + CZx²)` → 0).

The math reproduces the legacy `R!` exactly (same branch tolerance, same
special case), so transformation matrices match the pinned characterization
oracles to machine precision.
"""
function local_frame(x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real}, Ψ::Real; tol::Real=1e-6)
    v = SVector{3}(x2) - SVector{3}(x1)
    x = v / norm(v)
    CXx, CYx, CZx = x

    cΨ = cos(Ψ)
    sΨ = sin(Ψ)

    # cross(x, globalY) = (−CZx, 0, CXx); small norm ⇒ member parallel to global Y
    if sqrt(CZx^2 + CXx^2) < tol
        return SMatrix{3,3}(
            0.0, -CYx * cΨ, CYx * sΨ,   # column 1
            CYx, 0.0, 0.0,               # column 2
            0.0, sΨ, cΨ,                 # column 3
        )
    else
        d = sqrt(CXx^2 + CZx^2)
        b1 = (-CXx * CYx * cΨ - CZx * sΨ) / d
        b2 = d * cΨ
        b3 = (-CYx * CZx * cΨ + CXx * sΨ) / d
        c1 = (CXx * CYx * sΨ - CZx * cΨ) / d
        c2 = -d * sΨ
        c3 = (CYx * CZx * sΨ + CXx * cΨ) / d

        return SMatrix{3,3}(
            CXx, b1, c1,   # column 1
            CYx, b2, c2,   # column 2
            CZx, b3, c3,   # column 3
        )
    end
end

"""
    element_length(x1, x2) -> T

Distance between the element's end positions [length].
"""
element_length(x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real}) =
    norm(SVector{3}(x2) - SVector{3}(x1))

"""
    transform_to_global(k_local::SMatrix{12,12}, Λ::SMatrix{3,3}) -> SMatrix{12,12}

Rotate an element stiffness matrix from local to global coordinates:
`K = Rᵀ k R` where `R = blockdiag(Λ, Λ, Λ, Λ)`.

Computed blockwise over the sixteen 3×3 sub-blocks (`K[a,b] = Λᵀ k[a,b] Λ`)
— algebraically identical to the dense triple product but cheaper, and pure
SMatrix arithmetic so it is stack-allocated and AD-transparent.
"""
function transform_to_global(k::SMatrix{12,12,T}, Λ::SMatrix{3,3,TΛ}) where {T,TΛ}
    TK = promote_type(T, TΛ)
    K = _mutable12(TK)
    @inbounds for a in 0:3, c in 0:3
        ia = SVector(3a + 1, 3a + 2, 3a + 3)
        ic = SVector(3c + 1, 3c + 2, 3c + 3)
        K[ia, ic] = Λ' * k[ia, ic] * Λ
    end
    return SMatrix{12,12,TK}(K)
end
