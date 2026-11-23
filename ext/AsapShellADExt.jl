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

# one tag type for every shell rule: the dual kernels compile once
struct ShellTag end
const _CHUNK = 4   # dual width: 4 balances compile (≈85 s) vs runtime (12 → 155 s compile, 2 → 2× slower)

@inline _pack(A, D, S, γt, X1, X2, X3, X4) =
    vcat(X1, X2, X3, X4, vec(A), vec(D), vec(S), SVector(γt))

@inline function _unpack35(p)
    X1 = p[SVector(1, 2, 3)]; X2 = p[SVector(4, 5, 6)]
    X3 = p[SVector(7, 8, 9)]; X4 = p[SVector(10, 11, 12)]
    A = SMatrix{3,3}(p[SVector(13, 14, 15, 16, 17, 18, 19, 20, 21)])
    D = SMatrix{3,3}(p[SVector(22, 23, 24, 25, 26, 27, 28, 29, 30)])
    S = SMatrix{2,2}(p[SVector(31, 32, 33, 34)])
    return A, D, S, p[35], X1, X2, X3, X4
end

_sv3(x) = SVector{3}(x)

function ChainRulesCore.rrule(::typeof(Asap._shell_K_plain), A, D, S, γt, X1, X2, X3, X4)
    K = Asap._shell_K_plain(A, D, S, γt, X1, X2, X3, X4)
    function shell_K_pullback(ΔK)
        Δ = SMatrix{24,24}(unthunk(ΔK))
        p0 = Vector(_pack(SMatrix{3,3}(A), SMatrix{3,3}(D), SMatrix{2,2}(S), γt,
            _sv3(X1), _sv3(X2), _sv3(X3), _sv3(X4)))
        h(p) = sum(Δ .* Asap._shell_K_plain(_unpack35(SVector{35}(p))...))
        cfg = ForwardDiff.GradientConfig(nothing, p0, ForwardDiff.Chunk{_CHUNK}(), ShellTag())
        g = ForwardDiff.gradient(h, p0, cfg, Val(false))
        return (NoTangent(), SMatrix{3,3}(g[13:21]), SMatrix{3,3}(g[22:30]),
            SMatrix{2,2}(g[31:34]), g[35],
            SVector{3}(g[1:3]), SVector{3}(g[4:6]), SVector{3}(g[7:9]), SVector{3}(g[10:12]))
    end
    return K, shell_K_pullback
end

function ChainRulesCore.rrule(::typeof(Asap._shell_fef_plain), A, q, X1, X2, X3, X4)
    f = Asap._shell_fef_plain(A, q, X1, X2, X3, X4)
    function shell_fef_pullback(Δf)
        Δ = SVector{24}(unthunk(Δf))
        p0 = Vector(vcat(_sv3(X1), _sv3(X2), _sv3(X3), _sv3(X4), vec(SMatrix{3,3}(A)), _sv3(q)))
        function h(p)
            v = SVector{24}(p)
            X = (v[SVector(1, 2, 3)], v[SVector(4, 5, 6)], v[SVector(7, 8, 9)], v[SVector(10, 11, 12)])
            Am = SMatrix{3,3}(v[SVector(13, 14, 15, 16, 17, 18, 19, 20, 21)])
            return dot(Δ, Asap._shell_fef_plain(Am, v[SVector(22, 23, 24)], X...))
        end
        cfg = ForwardDiff.GradientConfig(nothing, p0, ForwardDiff.Chunk{_CHUNK}(), ShellTag())
        g = ForwardDiff.gradient(h, p0, cfg, Val(false))
        return (NoTangent(), SMatrix{3,3}(g[13:21]), SVector{3}(g[22:24]),
            SVector{3}(g[1:3]), SVector{3}(g[4:6]), SVector{3}(g[7:9]), SVector{3}(g[10:12]))
    end
    return f, shell_fef_pullback
end

end # module
