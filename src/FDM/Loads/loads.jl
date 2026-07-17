"""
    FDMload(points::Vector{FDMnode}, i::Int64, force::Vector{<:Real})
    FDMload(point::FDMnode, force::Vector{<:Real})

An external force applied to a FDM node: `load = [Px, Py, Pz]`
"""
mutable struct FDMload
    point::FDMnode # point at which load is applied
    force::Vector{Float64} # force vector

    function FDMload(points::Vector{FDMnode}, i::Int64, force::Vector{<:Real})
        @assert length(force) == 3 "Force vector should be length 3"

        load = new(points[i], Float64.(force))
        return load
    end

    function FDMload(point::FDMnode, force::Vector{<:Real})
        @assert length(force) == 3 "Force vector should be length 3"

        load = new(point, Float64.(force))
        return load
    end
end