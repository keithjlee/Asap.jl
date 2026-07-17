"""
    FDMelement{T}

An element in a Force Density Method network: it connects two nodes with a
force density `q = N/L` [force/length]. The member force follows the
geometry: `force(el) = q · length(el)`, always evaluated at the CURRENT node
positions.

`q` is mutable design data — change it freely and re-`solve!`; assembly
reads it fresh every time (see [`update_q!`](@ref) for the batch version).

# Fields
- `pStart`, `pEnd::FDMnode{T}`: the connected nodes
- `q::T`: force density [force/length]
- `id::Symbol`: identifier (optional)
- `iStart`, `iEnd::Int`: node indices (internal, assigned by `process!`)
- `index::Int`: element index (internal)

# Constructors
    FDMelement(pointStart, pointEnd, q, id = :element)
    FDMelement(points, iStart::Int, iEnd::Int, q, id = :element)
    FDMelement(points, indices::Vector{Int}, q, id = :element)
"""
mutable struct FDMelement{T<:Real}
    pStart::FDMnode{T}
    pEnd::FDMnode{T}
    q::T
    id::Symbol
    iStart::Int
    iEnd::Int
    index::Int

    function FDMelement(pointStart::FDMnode{T}, pointEnd::FDMnode{T}, q::Real,
        id::Symbol=:element) where {T}
        return new{T}(pointStart, pointEnd, T(q), id, 0, 0, 0)
    end
end

FDMelement(points::Vector{<:FDMnode}, iStart::Int, iEnd::Int, q::Real, id::Symbol=:element) =
    FDMelement(points[iStart], points[iEnd], q, id)
FDMelement(points::Vector{<:FDMnode}, indices::Vector{Int}, q::Real, id::Symbol=:element) =
    FDMelement(points[indices[1]], points[indices[2]], q, id)

include("utilities.jl")
