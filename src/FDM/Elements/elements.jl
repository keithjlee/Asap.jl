"""
    FDMelement(points::Vector{FDMnode}, iStart::Int64, iEnd::Int64, q::Real, id = :element)
    FDMelement(points::Vector{FDMnode}, indices::Vector{Int64}, q::Real, id = :element)
    FDMelement(pointStart::FDMnode, pointEnd::FDMnode, q::Real, id = :element)

An element in a Force Density Method network that connects two nodes with force density `q` and an optional identifier `id::Symbol`
"""
mutable struct FDMelement
    pStart::FDMnode #start point
    pEnd::FDMnode #end point
    q::Float64 #force density
    id::Symbol
    iStart::Int64 #index of start point in vector of points
    iEnd::Int64 #index of end point in vector of points
    elementID::Int64

    function FDMelement(points::Vector{FDMnode}, iStart::Int64, iEnd::Int64, q::Real, id = :element)
        element = new(points[iStart], points[iEnd], q, id)
        return element
    end

    function FDMelement(points::Vector{FDMnode}, indices::Vector{Int64}, q::Real, id = :element)
        iStart, iEnd = indices
        element = new(points[iStart], points[iEnd], q, id)
        return element
    end

    function FDMelement(pointStart::FDMnode, pointEnd::FDMnode, q::Real, id = :element)
        element = new(pointStart, pointEnd, q, id)
        return element
    end

    
end

include("utilities.jl")