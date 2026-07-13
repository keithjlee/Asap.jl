#=
Internal force (and displacement) recovery along elements — equilibrium
method.

The recovered section forces come from statics: starting from the EXACT
local end actions (f = k·u + q̃, which the finite element solution gets
exactly right at the nodes for Euler-Bernoulli members with consistent
fixed-end forces), the applied loads are integrated analytically along the
member. The result is exact at every station, for any end condition —
releases and semi-rigid springs need no special cases because they are
already baked into the end actions.

This replaces the per-(load type × release) analytic formula catalog that
lived in AsapToolkit's ForceAnalysis (~600 lines, incomplete, and with a
known axis-naming trap: its `.My` was the moment about local z). Here the
names are axis-correct: `Mz` pairs with `Vy` (bending in the local x–y
plane), `My` with `Vz`.

Sign conventions (local axes; x along the member):
- `N`: axial force, tension positive
- `Vy`, `Vz`: transverse shear — positive when the left piece pushes the
  right piece in the positive local axis direction
- `Mz`: bending moment about local z, sagging-positive for +y loads;
  satisfies dMz/dx = Vy
- `My`: bending moment about local y; satisfies dMy/dx = −Vz (the classic
  x–z plane sign flip)
- `Mx`: torsion (St. Venant twisting moment)
=#

"""
    LoadTrace{T}

The merged LOCAL loading of one prismatic member (or segment of a
super-element): everything needed to evaluate internal forces at any
station by closed form.

# Fields
- `L::T`: member length
- `xb::Vector{T}`: breakpoint stations `0 = xb[1] < … < xb[end] = L`
  (union of all distributed-load breakpoints)
- `w::Matrix{T}`: `3 × n` local distributed intensity components
  (wx, wy, wz) at each breakpoint, piecewise linear between [force/length]
- `W0::Matrix{T}`: cumulative `∫₀^xb w ds` per component (zeroth moment)
- `W1::Matrix{T}`: cumulative `∫₀^xb w·s ds` per component (first moment) —
  together these make `∫₀ˣ w(s)(x−s) ds = x·W0(x) − W1(x)` a closed form
- `pstation::Vector{T}`: point-action stations
- `pforce::Matrix{T}`: `3 × np` local point force components
- `pmoment::Matrix{T}`: `3 × np` local point moment components
"""
struct LoadTrace{T}
    L::T
    xb::Vector{T}
    w::Matrix{T}
    W0::Matrix{T}
    W1::Matrix{T}
    pstation::Vector{T}
    pforce::Matrix{T}
    pmoment::Matrix{T}
end

# piecewise-linear intensity component c at station x (exact)
function _trace_w(trace::LoadTrace{T}, c::Int, x::Real) where {T}
    xb = trace.xb
    k = clamp(searchsortedlast(xb, x), 1, length(xb) - 1)
    h = xb[k+1] - xb[k]
    h <= 0 && return trace.w[c, k+1]         # zero-width jump interval: right limit
    frac = (x - xb[k]) / h
    return trace.w[c, k] + (trace.w[c, k+1] - trace.w[c, k]) * frac
end

# cumulative ∫₀ˣ w_c ds and ∫₀ˣ w_c·s ds (exact: quadratic/cubic tails)
function _trace_integrals(trace::LoadTrace{T}, c::Int, x::Real) where {T}
    xb = trace.xb
    k = clamp(searchsortedlast(xb, x), 1, length(xb) - 1)
    x0 = xb[k]
    h = x - x0
    h <= 0 && return trace.W0[c, k], trace.W1[c, k]    # at a breakpoint (incl. jumps)
    w0 = trace.w[c, k]
    wx = _trace_w(trace, c, x)
    seg0 = h * (w0 + wx) / 2                                       # ∫ tail
    # ∫ s·w ds over [x0, x] with linear w: trapezoid on s·w is not exact —
    # use the exact cubic: ∫ s·(a + b·s) ds with a,b from the segment
    b = (wx - w0) / h
    a = w0 - b * x0
    seg1 = a * (x^2 - x0^2) / 2 + b * (x^3 - x0^3) / 3
    return trace.W0[c, k] + seg0, trace.W1[c, k] + seg1
end

"""
    ElementForceState{T}

Everything needed to evaluate a member's internal forces and local
displacements at any fraction of its length, in closed form and without
allocation: the local end actions and displacements at the start node, the
member's rigidities, and its merged [`LoadTrace`](@ref).

Build with [`internal_forces`](@ref); evaluate with [`axial_force`](@ref),
[`shear_y`](@ref)/[`shear_z`](@ref), [`moment_y`](@ref)/[`moment_z`](@ref),
[`torsion`](@ref), and [`local_displacements`](@ref); sample densely with
[`InternalForces`](@ref).
"""
struct ElementForceState{T}
    L::T
    f0::SVector{6,T}      # local end actions at the start node (on the element)
    u0::SVector{6,T}      # local displacements/rotations at the start node
    EA::T
    EIx::T
    EIy::T
    GJ::T
    trace::LoadTrace{T}
end

# ── scalar evaluators (t is a fraction of the member length) ────────────────

"""
    axial_force(state::ElementForceState, t) -> T

Internal axial force at fraction `t` [force], tension-positive.
"""
function axial_force(s::ElementForceState, t::Real)
    x = t * s.L
    W0x, _ = _trace_integrals(s.trace, 1, x)
    P = _point_sum(s, 1, x)
    return -(s.f0[1] + W0x + P)
end

"""
    shear_y(state::ElementForceState, t) -> T

Internal shear in the local y direction at fraction `t` [force].
Satisfies `dMz/dx = Vy`.
"""
function shear_y(s::ElementForceState, t::Real)
    x = t * s.L
    W0y, _ = _trace_integrals(s.trace, 2, x)
    return s.f0[2] + W0y + _point_sum(s, 2, x)
end

"""
    shear_z(state::ElementForceState, t) -> T

Internal shear in the local z direction at fraction `t` [force].
"""
function shear_z(s::ElementForceState, t::Real)
    x = t * s.L
    W0z, _ = _trace_integrals(s.trace, 3, x)
    return s.f0[3] + W0z + _point_sum(s, 3, x)
end

"""
    torsion(state::ElementForceState, t) -> T

Internal torsional (twisting) moment at fraction `t` [force·length].
"""
function torsion(s::ElementForceState, t::Real)
    x = t * s.L
    return -(s.f0[4] + _point_moment_sum(s, 1, x))
end

"""
    moment_z(state::ElementForceState, t) -> T

Internal bending moment about local z at fraction `t` [force·length] —
the moment paired with `Vy` (bending in the local x–y plane),
sagging-positive for +y loading.
"""
function moment_z(s::ElementForceState, t::Real)
    x = t * s.L
    W0y, W1y = _trace_integrals(s.trace, 2, x)
    lever = x * W0y - W1y                       # ∫ w_y(s)(x−s) ds
    return x * s.f0[2] - s.f0[6] + lever + _point_lever_sum(s, 2, x) - _point_moment_sum(s, 3, x)
end

"""
    moment_y(state::ElementForceState, t) -> T

Internal bending moment about local y at fraction `t` [force·length] —
paired with `Vz` (bending in the local x–z plane); `dMy/dx = −Vz`.
"""
function moment_y(s::ElementForceState, t::Real)
    x = t * s.L
    W0z, W1z = _trace_integrals(s.trace, 3, x)
    lever = x * W0z - W1z
    return -(x * s.f0[3] + s.f0[5] + lever + _point_lever_sum(s, 3, x) + _point_moment_sum(s, 2, x))
end

# Σ point-force component c applied at stations ≤ x
function _point_sum(s::ElementForceState{T}, c::Int, x::Real) where {T}
    acc = zero(T)
    @inbounds for (i, a) in enumerate(s.trace.pstation)
        a <= x && (acc += s.trace.pforce[c, i])
    end
    return acc
end

# Σ point-force component c times lever arm (x − a)
function _point_lever_sum(s::ElementForceState{T}, c::Int, x::Real) where {T}
    acc = zero(T)
    @inbounds for (i, a) in enumerate(s.trace.pstation)
        a <= x && (acc += s.trace.pforce[c, i] * (x - a))
    end
    return acc
end

# Σ point-moment component c at stations ≤ x
function _point_moment_sum(s::ElementForceState{T}, c::Int, x::Real) where {T}
    acc = zero(T)
    @inbounds for (i, a) in enumerate(s.trace.pstation)
        a <= x && (acc += s.trace.pmoment[c, i])
    end
    return acc
end

# ── displacement recovery ────────────────────────────────────────────────────

"""
    local_displacements(state::ElementForceState, t) -> SVector{3}

Local displacements `(u, v, w)` at fraction `t`: axial extension plus
transverse deflections in the local y and z directions [length].

Recovered exactly by double integration of the exact internal-force
fields — `v(x) = v₀ + ∫θz`, `θz(x) = θz₀ + ∫Mz/EIx` (and the x–z analogue
with the sign flip `θy = −w′`), axial `u(x) = u₀ + ∫N/EA`. The integrands
are piecewise polynomials of degree ≤ 3, so per-interval 3-point Gauss
evaluates the integrals exactly.
"""
function local_displacements(s::ElementForceState{T}, t::Real) where {T}
    x = t * s.L
    u = s.u0[1] + _gauss_integral(ξ -> axial_force(s, ξ / s.L) / s.EA, s, x)
    θz = ξ -> s.u0[6] + _gauss_integral(η -> moment_z(s, η / s.L) / s.EIx, s, ξ)
    v = s.u0[2] + _gauss_integral(θz, s, x)
    θy = ξ -> s.u0[5] + _gauss_integral(η -> moment_y(s, η / s.L) / s.EIy, s, ξ)
    w = s.u0[3] - _gauss_integral(θy, s, x)
    return SVector(u, v, w)
end

# ∫₀ˣ f dξ with breakpoint- and point-station-aware panels, 3-pt Gauss each
function _gauss_integral(f, s::ElementForceState{T}, x::Real) where {T}
    stations = _integration_stations(s, x)
    acc = zero(T)
    @inbounds for k in 1:length(stations)-1
        a, b = stations[k], stations[k+1]
        h = (b - a) / 2
        m = (a + b) / 2
        for g in 1:3
            acc += GAUSS3_W[g] * h * f(m + h * GAUSS3_X[g])
        end
    end
    return acc
end

function _integration_stations(s::ElementForceState{T}, x::Real) where {T}
    st = T[0]
    for b in s.trace.xb
        0 < b < x && push!(st, b)
    end
    for a in s.trace.pstation
        0 < a < x && push!(st, a)
    end
    push!(st, x)
    return sort!(unique!(st))
end

# ── construction ─────────────────────────────────────────────────────────────

"""
    internal_forces(model, el) -> ElementForceState
    internal_forces(model, el::VariableElement) -> Vector{ElementForceState}

Build the closed-form internal-force state of an element from a solved
model (requires `solve!` to have run). For a `VariableElement`, one state
per segment; use [`locate_segment`](@ref) or the fraction-based evaluators
below, which handle the mapping for you:

    axial_force(model, el, t), shear_y(model, el, t), moment_z(model, el, t), …
"""
function internal_forces(model::Model{T}, el::Union{FrameElement,TrussElement};
    results::Union{Nothing,LinearResults{T}}=nothing,
    factors::Union{Nothing,Dict{Symbol,T}}=nothing) where {T}
    res = results === nothing ? model.results::LinearResults{T} : results
    L = Base.length(el)
    Λ = local_frame(el)
    f = element_forces(res, el)
    u_el = displacement(res, el.nodeStart)
    u0 = vcat(Λ * SVector{3,T}(u_el[1], u_el[2], u_el[3]),
        Λ * SVector{3,T}(u_el[4], u_el[5], u_el[6]))
    sec = el.section
    trace = build_trace(model, el, Λ, L; factors=factors)
    return ElementForceState{T}(L, SVector{6,T}(f[1:6]), SVector{6,T}(u0),
        EA(sec), EIx(sec), EIy(sec), GJ(sec), trace)
end

function internal_forces(model::Model{T}, el::VariableElement;
    results::Union{Nothing,LinearResults{T}}=nothing,
    factors::Union{Nothing,Dict{Symbol,T}}=nothing) where {T}
    res = results === nothing ? model.results::LinearResults{T} : results
    Λ = local_frame(el)
    u = res.u
    g = element_global_dofs(el)

    return map(1:n_segments(el)) do s
        xa, xb_ = segment_positions(el, s)
        Ls = element_length(xa, xb_)
        slots = segment_slots(el, s)
        f = element_forces(res, el, s)
        us = SVector{6,T}(ntuple(i -> u[g[slots[i]]], 6))
        u0 = vcat(Λ * SVector{3,T}(us[1], us[2], us[3]),
            Λ * SVector{3,T}(us[4], us[5], us[6]))
        sec = el.sections[s]
        trace = build_trace(model, el, Λ, Ls; segment=s, factors=factors)
        ElementForceState{T}(Ls, SVector{6,T}(f[1:6]), SVector{6,T}(u0),
            EA(sec), EIx(sec), EIy(sec), GJ(sec), trace)
    end
end

"""
    build_trace(model, el, Λ, L; segment = nothing) -> LoadTrace

Merge every element load targeting `el` (restricted to `segment` for
super-elements) into one local [`LoadTrace`](@ref).
"""
function build_trace(model::Model{T}, el::AbstractElement, Λ::SMatrix{3,3}, L::T;
    segment::Union{Nothing,Int}=nothing,
    factors::Union{Nothing,Dict{Symbol,T}}=nothing) where {T}

    breaks = T[0, 1]
    dists = Tuple{Vector{T},Vector{T},SVector{3,T}}[]   # (t, w, local dir)
    pst = T[]
    pF = SVector{3,T}[]
    pM = SVector{3,T}[]

    for load in model.loads
        load isa ElementLoad || continue
        load.element === el || continue
        eff = segment === nothing ? load :
              _restrict_load(load, segment_fractions(el)[segment],
            segment_fractions(el)[segment+1], el)
        eff === nothing && continue
        λ = factors === nothing ? one(T) : get(factors, load.case, zero(T))
        iszero(λ) && continue
        _trace_add!(breaks, dists, pst, pF, pM, eff, Λ, el, segment, λ)
    end

    sort!(unique!(breaks))
    # DOUBLE every breakpoint: partial-span loads with nonzero boundary
    # intensity are JUMPS in w(x); a zero-width interval between duplicated
    # stations carries the jump exactly (left limit at the first copy, right
    # limit at the second), so cumulative integrals and evaluators stay exact
    doubled = T[]
    for t in breaks
        push!(doubled, t, t)
    end
    nb = Base.length(doubled)
    xb = doubled .* L
    w = zeros(T, 3, nb)
    for (td, wd, dloc) in dists
        for k in 1:nb
            side = isodd(k) ? :below : :above
            wt = _pw_linear_side(td, wd, doubled[k], side)
            for c in 1:3
                w[c, k] += wt * dloc[c]
            end
        end
    end

    # cumulative integrals at breakpoints (exact trapezoids of linear pieces)
    W0 = zeros(T, 3, nb)
    W1 = zeros(T, 3, nb)
    for k in 2:nb
        h = xb[k] - xb[k-1]
        for c in 1:3
            wa, wb = w[c, k-1], w[c, k]
            W0[c, k] = W0[c, k-1] + h * (wa + wb) / 2
            b = h > 0 ? (wb - wa) / h : zero(T)
            a = wa - b * xb[k-1]
            W1[c, k] = W1[c, k-1] + a * (xb[k]^2 - xb[k-1]^2) / 2 +
                       b * (xb[k]^3 - xb[k-1]^3) / 3
        end
    end

    np = Base.length(pst)
    pforce = zeros(T, 3, np)
    pmoment = zeros(T, 3, np)
    for i in 1:np
        pforce[:, i] = pF[i]
        pmoment[:, i] = pM[i]
    end

    return LoadTrace{T}(L, xb, w, W0, W1, pst .* L, pforce, pmoment)
end

# piecewise-linear load evaluation with one-sided limits at the support
# boundary: :below takes the limit from smaller t (zero AT the lower
# boundary), :above from larger t (zero AT the upper boundary)
function _pw_linear_side(td, wd, t, side::Symbol)
    z = zero(eltype(wd))
    if side === :below
        (t <= first(td) || t > last(td)) && return z
    else
        (t < first(td) || t >= last(td)) && return z
    end
    k = clamp(searchsortedlast(td, t), 1, Base.length(td) - 1)
    td[k+1] == td[k] && return wd[k]
    return wd[k] + (wd[k+1] - wd[k]) * (t - td[k]) / (td[k+1] - td[k])
end

function _trace_add!(breaks, dists, pst, pF, pM, load::DistributedLoad{T}, Λ, el, segment, λ) where {T}
    d = _local_direction(load.direction, load.coords, Λ)
    append!(breaks, load.t)
    push!(dists, (load.t, λ .* load.w, d))
    return
end

function _trace_add!(breaks, dists, pst, pF, pM, load::PointLoad{T}, Λ, el, segment, λ) where {T}
    p = load.coords === :local ? load.value : Λ * load.value
    push!(pst, load.position)
    push!(pF, λ * p)
    push!(pM, zero(SVector{3,T}))
    return
end

function _trace_add!(breaks, dists, pst, pF, pM, load::PointMoment{T}, Λ, el, segment, λ) where {T}
    m = load.coords === :local ? load.value : Λ * load.value
    push!(pst, load.position)
    push!(pF, zero(SVector{3,T}))
    push!(pM, λ * m)
    return
end

function _trace_add!(breaks, dists, pst, pF, pM, load::SelfWeight{T}, Λ, el, segment, λ) where {T}
    gmag = norm(load.g)
    iszero(gmag) && return
    sec = segment === nothing ? el.section : el.sections[segment]
    wmag = λ * ρA(sec) * gmag * load.factor
    d = Λ * (load.g ./ gmag)
    push!(dists, (T[0, 1], [wmag, wmag], d))
    return
end

# ── fraction-based convenience (uniform across primitive & super elements) ──

for f in (:axial_force, :shear_y, :shear_z, :torsion, :moment_y, :moment_z)
    @eval begin
        $f(model::Model, el::Union{FrameElement,TrussElement}, t::Real) =
            $f(internal_forces(model, el), t)
        function $f(model::Model, el::VariableElement, t::Real)
            s, τ = locate_segment(el, t)
            return $f(internal_forces(model, el)[s], τ)
        end
    end
end

# ── dense sampling for plotting ──────────────────────────────────────────────

"""
    InternalForces{T}

Densely sampled internal-force diagrams of one member, for plotting.
Stations include every load breakpoint and BOTH sides of each point action
(offset by ~√eps relative), so shear/moment discontinuities render as true
jumps instead of aliased slopes.

# Fields
- `x::Vector{T}`: stations along the member [length]
- `N`: axial force (tension +)
- `Vy`, `Mz`: shear and moment of the local x–y bending plane
- `Vz`, `My`: shear and moment of the local x–z bending plane
- `Mx`: torsion

(NOTE for AsapToolkit migrants: names are axis-correct here — the legacy
Toolkit `.My` corresponds to `Mz`, `.Mz` to `My`, `.P` to `N`.)
"""
struct InternalForces{T}
    x::Vector{T}
    N::Vector{T}
    Vy::Vector{T}
    Mz::Vector{T}
    Vz::Vector{T}
    My::Vector{T}
    Mx::Vector{T}
end

"""
    InternalForces(model, el; resolution = 20) -> InternalForces
    InternalForces(state::ElementForceState; resolution = 20)

Sample a member's internal forces at `resolution` evenly spaced stations
plus all load breakpoints and both sides of point actions. For a
`VariableElement`, segments are sampled in turn with globally increasing
stations.
"""
function InternalForces(state::ElementForceState{T}; resolution::Int=20) where {T}
    ts = _sample_fractions(state, resolution)
    return InternalForces{T}(
        ts .* state.L,
        [axial_force(state, t) for t in ts],
        [shear_y(state, t) for t in ts],
        [moment_z(state, t) for t in ts],
        [shear_z(state, t) for t in ts],
        [moment_y(state, t) for t in ts],
        [torsion(state, t) for t in ts],
    )
end

InternalForces(model::Model{T}, el::Union{FrameElement,TrussElement}; resolution::Int=20,
    results::Union{Nothing,LinearResults{T}}=nothing,
    factors::Union{Nothing,Dict{Symbol,T}}=nothing) where {T} =
    InternalForces(internal_forces(model, el; results=results, factors=factors); resolution=resolution)

function InternalForces(model::Model{T}, el::VariableElement; resolution::Int=20,
    results::Union{Nothing,LinearResults{T}}=nothing,
    factors::Union{Nothing,Dict{Symbol,T}}=nothing) where {T}
    states = internal_forces(model, el; results=results, factors=factors)
    ts = segment_fractions(el)
    L = Base.length(el)
    per = max(2, ceil(Int, resolution / n_segments(el)))
    xs = T[]
    cols = [T[] for _ in 1:6]
    for (s, st) in enumerate(states)
        for τ in _sample_fractions(st, per)
            push!(xs, (ts[s] + τ * (ts[s+1] - ts[s])) * L)
            push!(cols[1], axial_force(st, τ))
            push!(cols[2], shear_y(st, τ))
            push!(cols[3], moment_z(st, τ))
            push!(cols[4], shear_z(st, τ))
            push!(cols[5], moment_y(st, τ))
            push!(cols[6], torsion(st, τ))
        end
    end
    return InternalForces{T}(xs, cols...)
end

function _sample_fractions(state::ElementForceState{T}, resolution::Int) where {T}
    ts = collect(range(zero(T), one(T); length=max(resolution, 2)))
    ε = sqrt(eps(T))
    for b in state.trace.xb
        τ = b / state.L
        0 < τ < 1 && push!(ts, τ)
    end
    for a in state.trace.pstation
        τ = a / state.L
        0 < τ < 1 && append!(ts, (τ - ε, τ + ε))
    end
    return sort!(unique!(clamp.(ts, zero(T), one(T))))
end
