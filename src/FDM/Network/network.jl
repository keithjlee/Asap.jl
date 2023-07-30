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

    function Network(nodes::Vector{FDMnode}, elements::Vector{FDMelement}, loads::Vector{FDMload})

        network = new(nodes, elements, loads)

        network.processed = false
        return network
    end
end
