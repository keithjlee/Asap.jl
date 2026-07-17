"""
    Network{T}

A Force Density Method network: nodes, elements, and loads — pure
definition data, mirroring the frame core's three-layer separation. The
analysis structure ([`NetworkCache`](@ref), built by `process!`) and the
per-solve [`NetworkResults`](@ref) live in their own layers:

```julia
network = Network(nodes, elements, loads)
solve!(network)                        # geometry updates in place
member_force(network.results, el)
reaction(network.results, anchor)
```

Force densities, anchor positions, and loads are read FRESH at every
`solve!` — mutate them freely and re-solve at the cost of a numeric-only
refactorization on the frozen pattern. Adding/removing nodes, elements, or
loads, or changing fixity, requires `solve!(network; reprocess = true)`.

# Fields
- `nodes::Vector{FDMnode{T}}`
- `elements::Vector{FDMelement{T}}`
- `loads::Vector{FDMload{T}}`
- `cache::Union{Nothing,NetworkCache{T}}`: per-topology analysis structure
- `results::Union{Nothing,NetworkResults{T}}`: forces/reactions of the last solve
"""
mutable struct Network{T<:Real}
    nodes::Vector{FDMnode{T}}
    elements::Vector{FDMelement{T}}
    loads::Vector{FDMload{T}}
    cache::Any        # Union{Nothing,NetworkCache{T}} (cache type defined after Network)
    results::Any      # Union{Nothing,NetworkResults{T}}

    function Network(nodes::Vector{FDMnode{T}}, elements::Vector{FDMelement{T}},
        loads::Vector{FDMload{T}}) where {T}
        return new{T}(nodes, elements, loads, nothing, nothing)
    end
end

#narrowing convenience: accept loosely-typed containers (e.g. generators
#collecting into Vector{FDMelement}) and concretize on the node eltype
function Network(nodes::AbstractVector{<:FDMnode{T}}, elements::AbstractVector{<:FDMelement},
    loads::AbstractVector{<:FDMload}) where {T}
    return Network(collect(FDMnode{T}, nodes), collect(FDMelement{T}, elements),
        collect(FDMload{T}, loads))
end
