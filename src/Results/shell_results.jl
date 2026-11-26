#=
Shell results: element nodal forces, stress resultants, section cuts.
=#

# global node-major 6N displacement vector of an element
_element_u(res::LinearResults, el::AbstractElement) = res.u[element_global_dofs(el)]

# shell (any node count N): element nodal forces in LOCAL axes f = K_l·a_l + q̃
# (6N-vector), and the fixed-DOF entries accumulated into the reactions (global)
function _recover_forces!(reactions, isfixed, cache::AnalysisCache{T},
    el::ShellElement{T,N}, u::Vector{T}) where {T,N}
    X = map(n -> n.position, el.nodes)
    R = shell_frame(X...)[1]
    g = element_global_dofs(el)
    Ke = _shell_stiffness_d(el.section, X, el.drilling)
    Qg = _rotate_nodes(R', SVector{6N,T}(cache.q_local[el.index][1]))
    f_global = Ke * SVector{6N,T}(ntuple(i -> u[g[i]], Val(6N))) + Qg
    @inbounds for i in 1:6N
        isfixed[g[i]] && (reactions[g[i]] += f_global[i])
    end
    return Vector{T}(_rotate_nodes(R, f_global))
end

"""
    shell_resultants(model, el::ShellElement; at = (0.0, 0.0), results = model.results,
                     factors = nothing) -> NamedTuple

Stress resultants of a solved shell element at a point `at` (default: the
element centre), in the element's local axes ([`shell_frame`](@ref)):

- quadrilateral: natural coordinates `at = (ξ, η) ∈ [−1, 1]²`
- triangle: triangular coordinates `at = (r, s)` (node 1 at (0, 0), node 2 at
  (1, 0), node 3 at (0, 1); centroid (1/3, 1/3))

Fields:

- `N = (Nxx, Nyy, Nxy)` [force/length] — membrane forces, tension positive
- `M = (Mxx, Myy, Mxy)` [force·length/length] — bending/twisting moments
  per unit width; positive `Mxx` puts the `+e₃` face in tension along `e₁`
- `Q = (Qx, Qy)` [force/length] — transverse shear forces
- `frame::SMatrix{3,3}` — the local axes (rows `e₁, e₂, e₃`)

`results`/`factors` select which solution to read (e.g. a combination from
[`combine`](@ref) with its case factors), exactly as for
[`internal_forces`](@ref): `factors = Dict(case => factor)` weights the
element's own surface loads, which the membrane's condensed modes depend on.
Element centres are the superconvergent sampling points of the bilinear
fields; corner values are extrapolations — smooth them only deliberately.
"""
function shell_resultants(model::Model{T}, el::ShellElement{T,4};
    at = nothing, results::Union{Nothing,LinearResults{T}} = nothing,
    factors::Union{Nothing,Dict{Symbol,T}} = nothing) where {T}
    at = at === nothing ? (0.0, 0.0) : at
    res = results === nothing ? model.results::LinearResults{T} : results
    X = map(n -> n.position, el.nodes)
    R, xy, _ = shell_frame(X...)
    a = _rotate24(R, SVector{24,T}(_element_u(res, el)))
    fb = zero(SVector{4,T})
    for load in model.loads
        load isa ElementLoad && load.element === el || continue
        w = factors === nothing ? one(T) : get(factors, load.case, zero(T))
        iszero(w) && continue
        fb += w * shell_consistent_load(xy, R * _shell_traction(load, el))[2]
    end
    r = shell_resultants_local(el.section, xy, a, fb; at = at)
    return (N = r.N, M = r.M, Q = r.Q, frame = R)
end

function shell_resultants(model::Model{T}, el::ShellElement{T,3};
    at = nothing, results::Union{Nothing,LinearResults{T}} = nothing,
    factors::Union{Nothing,Dict{Symbol,T}} = nothing) where {T}
    at = at === nothing ? (1 / 3, 1 / 3) : at
    res = results === nothing ? model.results::LinearResults{T} : results
    X = map(n -> n.position, el.nodes)
    R, xy, _ = shell_frame(X...)
    a = _rotate18(R, SVector{18,T}(_element_u(res, el)))
    r = shell_resultants_local(el.section, xy, a; at = at,
        drilling = el.drilling == _ANDES_DRILL ? nothing : el.drilling)
    return (N = r.N, M = r.M, Q = r.Q, frame = R)
end

# global nodal forces (6 per node, nodes(el) order) of any supported element
function _global_nodal_forces(res::LinearResults{T}, el::ShellElement{T,N}) where {T,N}
    return _rotate_nodes(local_frame(el)', SVector{6N,T}(element_forces(res, el)))
end
function _global_nodal_forces(res::LinearResults{T}, el::Union{FrameElement,TrussElement}) where {T}
    return _rotate12(local_frame(el)', SVector{12,T}(element_forces(res, el)))
end
_global_nodal_forces(::LinearResults, el::AbstractElement) = throw(ArgumentError(
    "section_cut does not support $(nameof(typeof(el))) elements"))

"""
    section_cut(model, elements, cut_nodes; about = centroid(cut_nodes),
                axes = I, results = model.results) -> (F, M)

Resultant force `F` and moment `M` (about the point `about`) transmitted
ACROSS a cut into the free body formed by `elements`, through the nodes
`cut_nodes` that lie on the cut — the integrated "pier/spandrel" forces a
wall design check consumes.

Computed from the elements' nodal forces (`f = K·u + q̃`) at the cut nodes,
which makes it exact for equilibrium on any mesh: for a free body whose only
other actions are known loads, `F` and `M` balance them to machine
precision, independent of mesh refinement (no stress smoothing involved).
Frame and truss members crossing the cut (e.g. boundary columns embedded in
a wall) are included by listing them in `elements`.

- `elements`: elements on ONE side of the cut (any mix of `ShellElement`,
  `FrameElement`, `TrussElement`)
- `cut_nodes`: nodes on the cut line
- `about`: moment reference point (default: centroid of `cut_nodes`)
- `axes`: 3×3 matrix whose ROWS are the reporting directions (default:
  global); e.g. rows (wall length, wall thickness, vertical) give
  `F = (V₂, V₃, P)` and `M = (M₂, M₃, T)` for a pier
- `results`: any `LinearResults` (a case, or a `combine`d combination —
  element forces superpose)

Sign: `F`/`M` are the actions OF the rest of the structure ON `elements`.

# Example
```julia-repl
julia> below = filter(el -> centroid_z(el) < z_cut, wall_elements);

julia> cut = filter(n -> n.position[3] ≈ z_cut, wall_nodes);

julia> F, M = section_cut(model, below, cut)   # storey shear, axial, overturning
```
"""
function section_cut(model::Model{T}, elements::AbstractVector, cut_nodes::AbstractVector;
    about::Union{Nothing,AbstractVector} = nothing, axes::AbstractMatrix = SMatrix{3,3,T}(I),
    results::Union{Nothing,LinearResults{T}} = nothing) where {T}
    res = results === nothing ? model.results::LinearResults{T} : results
    isempty(cut_nodes) && throw(ArgumentError("section_cut needs at least one cut node"))
    cut = IdDict{Node{T},Bool}(n => true for n in cut_nodes)
    c = about === nothing ? sum(n -> n.position, cut_nodes) / length(cut_nodes) : SVector{3,T}(about)
    F = zero(SVector{3,T})
    M = zero(SVector{3,T})
    for el in elements
        f = _global_nodal_forces(res, el)
        for (k, n) in enumerate(nodes(el))
            haskey(cut, n) || continue
            fk = SVector{3,T}(f[6k-5], f[6k-4], f[6k-3])
            mk = SVector{3,T}(f[6k-2], f[6k-1], f[6k])
            F += fk
            M += mk + cross(n.position - c, fk)
        end
    end
    A = SMatrix{3,3,T}(axes)
    return (F = A * F, M = A * M)
end

"""
    shell_envelope(model, el::ShellElement, cr::CaseResults, combos; at = centre) -> NamedTuple

Component-wise extrema of the stress resultants of `el` at `at` over the load
combinations `combos` (the shell counterpart of [`envelope`](@ref)): fields
`Nmax`, `Nmin`, `Mmax`, `Mmin` (3-vectors, local `(xx, yy, xy)`), `Qmax`,
`Qmin` (2-vectors), each combination evaluated exactly by superposition
([`combine`](@ref)) with its case factors.
"""
function shell_envelope(model::Model{T}, el::ShellElement, cr::CaseResults{T},
    combos::Vector{<:LoadCombination}; at = nothing) where {T}
    isempty(combos) && throw(ArgumentError("need at least one combination to envelope"))
    rs = [shell_resultants(model, el; at = at, results = combine(cr, c), factors = factor_dict(c))
          for c in combos]
    pick(f, sym) = reduce((a, b) -> f.(a, b), getproperty.(rs, sym))
    return (Nmax = pick(max, :N), Nmin = pick(min, :N), Mmax = pick(max, :M),
            Mmin = pick(min, :M), Qmax = pick(max, :Q), Qmin = pick(min, :Q))
end
