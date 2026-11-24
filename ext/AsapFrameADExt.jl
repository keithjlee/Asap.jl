"""
    AsapFrameADExt

Reverse-mode rules for the frame kernels on plain data (`_frame_K_plain`,
`_fef_local_plain`, `_fef_global_plain`). Zygote tracing through a 12×12
static-matrix kernel records hundreds of scalar operations per element
(≈150 µs per stiffness pullback, ≈250 µs per element-load pullback, against
0.8 µs and 0.1 µs for the kernels themselves). These rules compute the
kernel's derivative with ForwardDiff instead — a few dual-number evaluations
of the SAME kernel — and contract it with the cotangent. Not a hand-derived
adjoint: the derivative is the kernel's own, in forward mode (the pattern of
`AsapShellADExt`).

Only Zygote (ChainRules) uses these; Enzyme and Mooncake differentiate the
frame kernels natively. Loads are data, never design variables: the load
argument gets no tangent. Active when ChainRulesCore and ForwardDiff are
both loaded (without ForwardDiff, reverse AD still traces the kernels).
"""
module AsapFrameADExt

using Asap, ChainRulesCore, ForwardDiff, StaticArrays, LinearAlgebra

struct FrameTag end
const _CHUNK = 10   # 19–20 inputs → 2 dual passes

_cotv(::Val{N}, Δ::AbstractArray) where {N} = SVector{N}(vec(Δ))
_cotv(v::Val, Δ::ChainRulesCore.Tangent) = _cotv(v, ChainRulesCore.backing(Δ))
_cotv(v::Val, Δ::NamedTuple) = _cotv(v, Δ.data)
_cotv(::Val{N}, Δ::Tuple) where {N} = SVector{N}(Δ)

const F64V{N} = SVector{N,Float64}

# ⟨Δ, f(p)⟩ gradient over the packed inputs, forward mode
function _vjp(h, p0::SVector{N,Float64}) where {N}
    v = Vector(p0)
    cfg = ForwardDiff.GradientConfig(nothing, v, ForwardDiff.Chunk{min(_CHUNK, N)}(), FrameTag())
    return SVector{N}(ForwardDiff.gradient(h, v, cfg, Val(false)))
end

# ── stiffness: inputs (r 4, k 8, x1 3, x2 3, rollangle 1) = 19 ──────────────
@inline _unpack_K(p) = (p[SVector(1, 2, 3, 4)], p[SVector(5, 6, 7, 8, 9, 10, 11, 12)],
    p[SVector(13, 14, 15)], p[SVector(16, 17, 18)], p[19])

function ChainRulesCore.rrule(::typeof(Asap._frame_K_plain), r::F64V{4}, k::F64V{8},
        x1::F64V{3}, x2::F64V{3}, rollangle::Float64)
    K = Asap._frame_K_plain(r, k, x1, x2, rollangle)
    function frame_K_pullback(ΔK)
        Δ = unthunk(ΔK)
        Δ isa ChainRulesCore.AbstractZero && return (NoTangent(), Δ, Δ, Δ, Δ, Δ)
        D = reshape(_cotv(Val(144), Δ), Size(12, 12))
        function h(p)
            rq, kq, a, b, θ = _unpack_K(SVector{19}(p))
            return sum(D .* Asap._frame_K_plain(rq, _hold_limits(kq, k), a, b, θ))
        end
        g = _vjp(h, vcat(r, k, x1, x2, SVector(rollangle)))
        gr, gk, g1, g2, gψ = _unpack_K(g)
        return (NoTangent(), gr, _finite(gk, k), g1, g2, gψ)
    end
    return K, frame_K_pullback
end

# ── element loads: inputs (r 3, k 8, x1 3, x2 3, rollangle 1) = 18 ──────────
@inline _unpack_F(p) = (p[SVector(1, 2, 3)], p[SVector(4, 5, 6, 7, 8, 9, 10, 11)],
    p[SVector(12, 13, 14)], p[SVector(15, 16, 17)], p[18])

for f in (:_fef_local_plain, :_fef_global_plain)
    @eval function ChainRulesCore.rrule(::typeof(Asap.$f), load, r::F64V{3}, k::F64V{8},
            x1::F64V{3}, x2::F64V{3}, rollangle::Float64)
        y = Asap.$f(load, r, k, x1, x2, rollangle)
        function fef_pullback(Δy)
            Δ = unthunk(Δy)
            Δ isa ChainRulesCore.AbstractZero && return (NoTangent(), NoTangent(), Δ, Δ, Δ, Δ, Δ)
            D = _cotv(Val(12), Δ)
            function h(p)
                rq, kq, a, b, θ = _unpack_F(SVector{18}(p))
                return dot(D, Asap.$f(load, rq, _hold_limits(kq, k), a, b, θ))
            end
            g = _vjp(h, vcat(r, k, x1, x2, SVector(rollangle)))
            gr, gk, g1, g2, gψ = _unpack_F(g)
            return (NoTangent(), NoTangent(), gr, _finite(gk, k), g1, g2, gψ)
        end
        return y, fef_pullback
    end
end

# Rigid (Inf) and released (0) springs are EXACT limits: the kernel takes them
# through `isinf`/`iszero` branches, so the plain trace gives them a zero
# derivative. Seeding a dual partial on one breaks that: an Inf value makes
# 0·Inf = NaN in every direction, and `iszero(::Dual)` also requires zero
# partials, so a seeded 0 would skip the release branch and divide by zero.
# The dual evaluation therefore holds limit springs CONSTANT, and their
# cotangent is zero, as in the plain trace.
_islimit(k) = !isfinite(k) || iszero(k)
_hold_limits(kq::SVector{N}, k::SVector{N}) where {N} =
    map((qi, ki) -> _islimit(ki) ? oftype(qi, ki) : qi, kq, k)
_finite(g::SVector{N}, k::SVector{N}) where {N} = map((gi, ki) -> _islimit(ki) ? zero(gi) : gi, g, k)

end # module
