abstract type AbstractGeo end

include("TrussGeometry.jl")
include("ModelGeometry.jl")
include("NetworkGeometry.jl")

"""
    Geo(model::Model) -> ModelGeo
    Geo(network::Network) -> NetworkGeo

Flatten a solved structure into plain plot-ready arrays: node positions,
displaced positions, element connectivity, and per-element force ranges with
their maxima. Dispatches on the structure type — a frame/truss [`Model`](@ref)
produces a [`ModelGeo`](@ref), an FDM [`Network`](@ref) a [`NetworkGeo`](@ref).
For an axial-only extraction of a truss model, construct a [`TrussGeo`](@ref)
directly.
"""
function Geo end

Geo(model::Model) = ModelGeo(model)
Geo(network::Network) = NetworkGeo(network)
