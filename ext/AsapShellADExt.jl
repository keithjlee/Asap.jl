"""
    AsapShellADExt

Reverse-mode rules for the flat-shell kernels. Zygote tracing through the
24×24 static-matrix kernels records thousands of scalar operations per
element (≈10 ms per stiffness pullback, minutes of compile); these rules
instead compute the kernel's local Jacobian with ForwardDiff — a handful of
dual-number evaluations of the SAME kernel — and contract it with the
cotangent. Not a hand-derived adjoint: the derivative is the kernel's own,
evaluated in forward mode.

Active when both ChainRulesCore and ForwardDiff are loaded (without
ForwardDiff, reverse AD still works through the generic trace, just slowly).
"""
module AsapShellADExt

using Asap, ChainRulesCore, ForwardDiff, StaticArrays, LinearAlgebra

# one tag type for every shell rule: the dual kernels compile once per width
struct ShellTag end

# Each pullback is the gradient of ⟨Δ, kernel(p)⟩ over the kernel's inputs,
# taken in SEPARATE ForwardDiff passes per input group: when only the
# rigidities (or the load) carry partials the geometry (shell frame,
# Jacobians, strain operators) stays plain Float64, which is most of the
# kernel. Only the position group pays for dual geometry.
# Dual width 4 for every group. Wider duals (one pass per group) run ~20% faster
# but compile far slower: first pullback 69 s vs 37 s on Julia 1.13, and 853 s vs
# 141 s on Julia 1.10, whose compiler struggles with wide duals through these
# unrolled 24×24 kernels.
const _CHUNK = 4

function _grad(h, p0::SVector{N,Float64}, ::Val{C}) where {N,C}
    v = Vector(p0)
    cfg = ForwardDiff.GradientConfig(nothing, v, ForwardDiff.Chunk{C}(), ShellTag())
    return SVector{N}(ForwardDiff.gradient(h, v, cfg, Val(false)))
end

@inline _m3(v, o) = SMatrix{3,3}(ntuple(k -> v[o+k], Val(9)))
@inline _x4(v) = (SVector(v[1], v[2], v[3]), SVector(v[4], v[5], v[6]),
    SVector(v[7], v[8], v[9]), SVector(v[10], v[11], v[12]))

_sv3(x) = SVector{3}(x)

# cotangent → flat static vector of length N, whatever the AD engine hands us:
# an array (Zygote), or a structural Tangent of a static array whose backing
# holds the `data` tuple (Mooncake's ChainRules interop)
_cot(::Val{N}, Δ::AbstractArray) where {N} = SVector{N}(vec(Δ))
_cot(v::Val, Δ::ChainRulesCore.Tangent) = _cot(v, ChainRulesCore.backing(Δ))   # may nest (Mooncake)
_cot(v::Val, Δ::NamedTuple) = _cot(v, Δ.data)
_cot(::Val{N}, Δ::Tuple) where {N} = SVector{N}(Δ)

function ChainRulesCore.rrule(::typeof(Asap._shell_K_plain), A, D, S, γt, X1, X2, X3, X4)
    K = Asap._shell_K_plain(A, D, S, γt, X1, X2, X3, X4)
    function shell_K_pullback(ΔK)
        Δ = reshape(_cot(Val(576), unthunk(ΔK)), Size(24, 24))
        As, Ds, Ss, γ = SMatrix{3,3,Float64}(A), SMatrix{3,3,Float64}(D), SMatrix{2,2,Float64}(S), Float64(γt)
        X = (_sv3(X1), _sv3(X2), _sv3(X3), _sv3(X4))
        # membrane (A, γt) and plate (D, S) at plain geometry
        gm = _grad(vcat(vec(As), SVector(γ)), Val(_CHUNK)) do p
            sum(Δ .* Asap._shell_K_plain(_m3(p, 0), Ds, Ss, p[10], X...))
        end
        gp = _grad(vcat(vec(Ds), vec(Ss)), Val(_CHUNK)) do p
            sum(Δ .* Asap._shell_K_plain(As, _m3(p, 0), SMatrix{2,2}(p[10], p[11], p[12], p[13]), γ, X...))
        end
        gx = _grad(vcat(X...), Val(_CHUNK)) do p
            sum(Δ .* Asap._shell_K_plain(As, Ds, Ss, γ, _x4(p)...))
        end
        x1, x2, x3, x4 = _x4(gx)
        return (NoTangent(), _m3(gm, 0), _m3(gp, 0), SMatrix{2,2}(gp[10], gp[11], gp[12], gp[13]), gm[10],
            x1, x2, x3, x4)
    end
    return K, shell_K_pullback
end

function ChainRulesCore.rrule(::typeof(Asap._shell_fef_plain), A, qg, ql, X1, X2, X3, X4)
    f = Asap._shell_fef_plain(A, qg, ql, X1, X2, X3, X4)
    function shell_fef_pullback(Δf)
        Δ = _cot(Val(24), unthunk(Δf))
        As, g0, l0 = SMatrix{3,3,Float64}(A), SVector{3,Float64}(qg), SVector{3,Float64}(ql)
        X = (_sv3(X1), _sv3(X2), _sv3(X3), _sv3(X4))
        # rigidity and traction at plain geometry; positions
        ga = _grad(vcat(vec(As), g0, l0), Val(_CHUNK)) do p
            dot(Δ, Asap._shell_fef_plain(_m3(p, 0), SVector(p[10], p[11], p[12]),
                SVector(p[13], p[14], p[15]), X...))
        end
        gx = _grad(vcat(X...), Val(_CHUNK)) do p
            dot(Δ, Asap._shell_fef_plain(As, g0, l0, _x4(p)...))
        end
        x1, x2, x3, x4 = _x4(gx)
        return (NoTangent(), _m3(ga, 0), SVector(ga[10], ga[11], ga[12]),
            SVector(ga[13], ga[14], ga[15]), x1, x2, x3, x4)
    end
    return f, shell_fef_pullback
end

# ShellSection → (A, D, S, γt): one width-4 dual pass over (E, G, ν, t) instead
# of reverse-tracing the four accessors (≈35 µs per element under Zygote)
_dense(Δ, y::Real) = Δ isa AbstractZero ? zero(y) : Δ
_dense(Δ, y::StaticArray) = Δ isa AbstractZero ? zero(y) : reshape(_cot(Val(length(y)), Δ), Size(y))

function ChainRulesCore.rrule(::typeof(Asap._shell_rigidities), s::ShellSection{Float64})
    y = Asap._shell_rigidities(s)
    function shell_rigidities_pullback(Δy)
        Δ = unthunk(Δy)
        Δ isa AbstractZero && return (NoTangent(), Δ)
        ΔA, ΔD, ΔS, Δγ = map(_dense, (unthunk(Δ[1]), unthunk(Δ[2]), unthunk(Δ[3]), unthunk(Δ[4])), y)
        m = s.material
        g = _grad(SVector(m.E, m.G, m.ν, s.t), Val(_CHUNK)) do p
            A, D, S, γ = Asap._shell_rigidities(ShellSection(Material(p[1], p[2], m.ρ, p[3]), p[4]))
            sum(ΔA .* A) + sum(ΔD .* D) + sum(ΔS .* S) + Δγ * γ
        end
        Δm = Tangent{typeof(m)}(E = g[1], G = g[2], ν = g[3])
        return (NoTangent(), Tangent{typeof(s)}(material = Δm, t = g[4]))
    end
    return y, shell_rigidities_pullback
end

end # module
