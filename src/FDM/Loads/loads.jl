"""
    FDMload(points::Vector{FDMnode}, i::Int64, force::Vector{<:Real})
    FDMload(point::FDMnode, force::Vector{<:Real})

An external force applied to a FDM node: `load = [Px, Py, Pz]`
"""
mutable struct FDMload
    point::FDMnode # point at which load is applied
    force::Vector{<:Real} # force vector

    function FDMload(points::Vector{FDMnode}, i::Int64, force::Vector{<:Real})
        @assert length(force) == 3 "Force vector should be length 3"

        load = new(points[i], force)
        return load
    end

    function FDMload(point::FDMnode, force::Vector{<:Real})
        @assert length(force) == 3 "Force vector should be length 3"

        load = new(point, force)
        return load
    end
end