"""
    FDMload{T}

An external force applied to an FDM node.

# Fields
- `point::FDMnode{T}`: the loaded node
- `force::SVector{3,T}`: the applied force [Px, Py, Pz] [force]

# Constructors
    FDMload(point, force)
    FDMload(points, i::Int, force)
"""
struct FDMload{T<:Real}
    point::FDMnode{T}
    force::SVector{3,T}

    function FDMload(point::FDMnode{T}, force::AbstractVector{<:Real}) where {T}
        @assert Base.length(force) == 3 "Force vector should be length 3"
        return new{T}(point, SVector{3,T}(force))
    end
end

FDMload(points::Vector{<:FDMnode}, i::Int, force::AbstractVector{<:Real}) =
    FDMload(points[i], force)
