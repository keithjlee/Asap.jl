#=
Load cases, combinations, and envelopes.

Linear analysis is linear: solve each load CASE once against a single
factorization (multi-RHS), then any weighted COMBINATION of cases is a
superposition — no re-solve, ever. Envelopes are extrema of combined
recovery over the combination set.
=#

"""
    LoadCombination{T}

A named, factored combination of load cases: `Σ factor · case`.

# Fields
- `name::Symbol`
- `factors::Vector{Pair{Symbol,T}}`: e.g. `[:dead => 1.2, :live => 1.6]`

# Examples
```julia-repl
julia> strength = LoadCombination(:LRFD1, [:dead => 1.2, :live => 1.6])

julia> service = LoadCombination(:service, [:dead => 1.0, :live => 1.0])
```
"""
struct LoadCombination{T<:Real}
    name::Symbol
    factors::Vector{Pair{Symbol,T}}
end

LoadCombination(name::Symbol, factors::Vector{<:Pair{Symbol,<:Real}}) =
    LoadCombination{Float64}(name, [Symbol(f.first) => Float64(f.second) for f in factors])

factor_dict(c::LoadCombination{T}) where {T} = Dict{Symbol,T}(c.factors)

"""
    CaseResults{T}

Results of a multi-case solve ([`solve_cases!`](@ref)): one
[`LinearResults`](@ref) per load case, all obtained from a single stiffness
assembly and a single factorization. Combine into any
[`LoadCombination`](@ref)'s results with [`combine`](@ref) — pure
superposition, no re-solve.

# Fields
- `cases::Vector{Symbol}`: the case tags, in solve order
- `results::Vector{LinearResults{T}}`: per-case results
- `F::Matrix{T}`: per-case full-space load vectors `P − Pf` (columns) —
  needed to compute combined compliance
"""
struct CaseResults{T}
    cases::Vector{Symbol}
    results::Vector{LinearResults{T}}
    F::Matrix{T}
end

Base.getindex(cr::CaseResults, case::Symbol) = cr.results[findfirst(==(case), cr.cases)]

"""
    load_cases(model) -> Vector{Symbol}

The distinct case tags present among the model's loads, in first-appearance
order.
"""
function load_cases(model::Model)
    cases = Symbol[]
    for load in model.loads
        load.case in cases || push!(cases, load.case)
    end
    return cases
end

"""
    solve_cases!(model; cases = load_cases(model)) -> CaseResults

Solve every load case against ONE stiffness assembly and ONE factorization:
the expensive work is done once, each case costs a load-vector assembly and
a pair of triangular back-substitutions. Combinations then come free via
[`combine`](@ref).
"""
function solve_cases!(model::Model{T}; cases::Vector{Symbol}=load_cases(model),
    reprocess::Bool=false) where {T}
    (model.cache === nothing || reprocess) && process!(model)
    cache = model.cache::AnalysisCache{T}

    assemble_K!(cache, model)
    cache.factorization = _factorize(cache.K)
    free = cache.partition.free

    results = LinearResults{T}[]
    F = zeros(T, cache.partition.n_global, Base.length(cases))
    for (j, case) in enumerate(cases)
        assemble_loads!(cache, model; case=case)
        F[:, j] = cache.P .- cache.Pf
        uf = cache.factorization \ F[free, j]
        u = zeros(T, cache.partition.n_global)
        u[free] = uf
        push!(results, _post_process(cache, model, u, dot(uf, F[free, j])))
    end

    # restore the all-cases load state so later solve!/recovery see everything
    assemble_loads!(cache, model)

    return CaseResults{T}(collect(cases), results, F)
end

"""
    combine(cr::CaseResults, combo::LoadCombination) -> LinearResults

Results of a factored combination by superposition: displacements,
reactions, and element end forces are linear in the loading, so the
combined results are exact — no additional solve. (Compliance, being
quadratic, is recomputed from the combined vectors.)

Cases named in the combination but absent from the solve contribute
nothing; a warning-free contract — validate combination spelling upstream
if needed.
"""
function combine(cr::CaseResults{T}, combo::LoadCombination) where {T}
    f = factor_dict(combo)
    λ = [get(f, c, zero(T)) for c in cr.cases]

    u = sum(λ[j] * cr.results[j].u for j in eachindex(λ))
    reactions = sum(λ[j] * cr.results[j].reactions for j in eachindex(λ))
    forces = [sum(λ[j] * cr.results[j].element_forces[i] for j in eachindex(λ))
              for i in eachindex(cr.results[1].element_forces)]
    Fc = cr.F * λ
    return LinearResults{T}(u, reactions, forces, dot(u, Fc))
end

"""
    Envelope{T}

Station-wise extrema of internal forces over a set of load combinations —
the quantity a design check actually consumes.

# Fields
- `x::Vector{T}`: stations along the member [length]
- `lo`, `hi`: `6 × n` matrices of minima/maxima, rows ordered
  (N, Vy, Mz, Vz, My, Mx)
- `combos::Vector{Symbol}`: names of the enveloped combinations
"""
struct Envelope{T}
    x::Vector{T}
    lo::Matrix{T}
    hi::Matrix{T}
    combos::Vector{Symbol}
end

"""
    envelope(model, el, cr::CaseResults, combos; resolution = 20) -> Envelope

Envelope a member's internal forces over the given combinations: each
combination's diagrams come from superposed case results (single
factorization behind all of it), and the envelope is their station-wise
extrema.
"""
function envelope(model::Model{T}, el::AbstractElement, cr::CaseResults{T},
    combos::Vector{<:LoadCombination}; resolution::Int=20) where {T}
    @assert !isempty(combos) "need at least one combination to envelope"

    # common stations for ALL combinations: derived with every case active
    # (factors = 1), so no combo's jump stations are missed even when some
    # combos zero out the loads that create them
    all_on = Dict{Symbol,T}(c => one(T) for c in cr.cases)
    ref = InternalForces(model, el; resolution=resolution,
        results=cr.results[1], factors=all_on)
    x = ref.x

    n = Base.length(x)
    lo = fill(T(Inf), 6, n)
    hi = fill(T(-Inf), 6, n)
    for combo in combos
        res = combine(cr, combo)
        f = _sample_at(model, el, x, res, factor_dict(combo))
        for (r, col) in enumerate(f)
            @inbounds for k in 1:n
                lo[r, k] = min(lo[r, k], col[k])
                hi[r, k] = max(hi[r, k], col[k])
            end
        end
    end
    return Envelope{T}(x, lo, hi, [c.name for c in combos])
end

# evaluate all six internal-force quantities at the given stations [length]
function _sample_at(model::Model{T}, el::Union{FrameElement,TrussElement},
    x::Vector{T}, res::LinearResults{T}, factors::Dict{Symbol,T}) where {T}
    st = internal_forces(model, el; results=res, factors=factors)
    ts = x ./ st.L
    return ([axial_force(st, t) for t in ts], [shear_y(st, t) for t in ts],
        [moment_z(st, t) for t in ts], [shear_z(st, t) for t in ts],
        [moment_y(st, t) for t in ts], [torsion(st, t) for t in ts])
end

function _sample_at(model::Model{T}, el::VariableElement,
    x::Vector{T}, res::LinearResults{T}, factors::Dict{Symbol,T}) where {T}
    states = internal_forces(model, el; results=res, factors=factors)
    L = Base.length(el)
    evals = ntuple(_ -> T[], 6)
    for xi in x
        s, τ = locate_segment(el, clamp(xi / L, zero(T), one(T)))
        st = states[s]
        push!(evals[1], axial_force(st, τ))
        push!(evals[2], shear_y(st, τ))
        push!(evals[3], moment_z(st, τ))
        push!(evals[4], shear_z(st, τ))
        push!(evals[5], moment_y(st, τ))
        push!(evals[6], torsion(st, τ))
    end
    return evals
end
