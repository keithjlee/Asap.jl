"""
FDMnode in fdm network. Includes positional data and fixed/free information
"""
mutable struct FDMnode
    x::Real
    y::Real
    z::Real
    dof::Bool # true = free; false = fixed
    id::Union{Symbol, Nothing}
    nodeID::Integer
    reaction::Vector{Float64}

    #empty constructor
    function FDMnode()
        return new()
    end

    # individual coordinate basis
    function FDMnode(x::Real, y::Real, z::Real, dof::Bool)
        return new(x, y, z, dof, nothing)
    end

    # using a vector to represent position
    function FDMnode(pos::Vector{<:Real}, dof::Bool)
        @assert length(pos) == 3 "pos should be length 3"
        
        return new(pos..., dof, nothing)
    end

    #with an id
    function FDMnode(x::Real, y::Real, z::Real, dof::Bool, id::Symbol)
        return new(x, y, z, dof, id)
    end

    function FDMnode(pos::Vector{<:Real}, dof::Bool, id::Symbol)
        @assert length(pos) == 3 "pos should be length 3"

        return new(pos..., dof, id)
    end
end

"""
FDMelement represents an edge in the network. Includes the global indices and positional information of the start and end nodes, as well as the force density.
"""
mutable struct FDMelement
    pStart::FDMnode #start point
    iStart::Int64 #index of start point in vector of points
    pEnd::FDMnode #end point
    iEnd::Int64 #index of end point in vector of points
    q::Float64 #force density
    id::Union{Symbol, Nothing}
    elementID::Integer

    function FDMelement(points::Vector{FDMnode}, iStart::Int64, iEnd::Int64, q::Real; id = nothing)
        element = new(points[iStart], iStart, points[iEnd], iEnd, Float64(q), id)
        return element
    end

    function FDMelement(points::Vector{FDMnode}, iStart::Int64, iEnd::Int64, q::Real, id::Symbol)
        element = new(points[iStart], iStart, points[iEnd], iEnd, Float64(q), id)
        return element
    end
end

"""
Load type is associated with a node (and its position in node vector), as well as a length 3 force vector (x,y,z)
"""
mutable struct FDMload
    point::FDMnode # point at which load is applied
    index::Int64 # position of point in vector of points
    force::Vector{<:Real} # force vector

    function FDMload(points::Vector{FDMnode}, i::Int64, force::Vector{<:Real})
        @assert length(force) == 3 "Force vector should be length 3"

        load = new(points[i], i, force)
        return load
    end
end

"""
An FDM network with all relevant information for analysis
"""
mutable struct Network
    nodes::Vector{FDMnode} #vector of nodes
    elements::Vector{FDMelement} #vector of elements
    loads::Vector{FDMload} #vector of loads
    q::Vector{<:Real} #vector of force densities
    Q::SparseMatrixCSC{Float64, Int64} #diagm(q)
    C::SparseMatrixCSC{Int64, Int64} #branch node matrix
    N::Vector{Int64} #fixed node indices
    F::Vector{Int64} #free node indices
    Cn::SparseMatrixCSC{Int64, Int64} #fixed branch node matrix
    Cf::SparseMatrixCSC{Int64, Int64} #free branch node matrix
    P::Matrix{<:Real}#load matrix
    Pn::Matrix{<:Real}#free node load matrix
    xyz::Union{Matrix{Int64}, Matrix{Float64}} #xyz matrix of initial positions
    processed::Bool

    function Network(nodes::Vector{FDMnode}, elements::Vector{FDMelement}, loads::Vector{FDMload}; copy = false)
        if copy
            network = new(deepcopy(nodes), deepcopy(elements), deepcopy(loads))
        else
            network = new(nodes, elements, loads)
        end
        network.processed = false
        return network
    end
end
