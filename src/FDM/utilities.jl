"""
Euclidean distance metric for FDMnode types
"""
function dist(a::FDMnode, b::FDMnode)
    return sqrt((a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2)
end

"""
Extract member forces
"""
function memberForces(network::Network)
    return norm.(eachrow(network.C * network.xyz)) .* network.q
end

"""
Extract member lengths
"""
function memberLengths(network::Network)
    return norm.(eachrow(network.C * network.xyz))
end

"""
Custom indexing based on IDs of structs
"""
function Base.getindex(nodes::Vector{FDMnode}, i::Symbol)
    return [node for node in nodes if node.id == i]
end

"""
Custom indexing for elements
"""
function Base.getindex(elements::Vector{FDMelement}, i::Symbol)
    return [element for element in elements if element.id == i]
end

"""
Custom indexing of networks
"""
function Base.getindex(network::Network, i::Symbol)
    nodes = network.nodes[i]
    elements = network.elements[i]

    return nodes, elements
end

"""
findall methods for querying
"""
function Base.findall(elements::Vector{FDMelement}, i::Symbol)
    return findall([x.id == i for x in elements])
end

function Base.findall(nodes::Vector{FDMnode}, i::Symbol)
    return findall([x.id == i for x in nodes])
end

function Base.findall(network::Network, i::Symbol)
    return findall(network.elements, i), findall(network.nodes, i)
end

"""
Repopulate network with new q values
"""
function qUpdate!(network::Network, q::Vector{<:Real})
    @assert length(network.elements) == length(q) "Number of elements and q must be equal"

    for (i, element) in enumerate(network.elements)
        element.q = q[i]
    end

    forceDensities!(network)
end

#single value
function qUpdate!(network::Network, q::Real)
    for element in network.elements
        element.q = q
    end

    forceDensities!(network)
end

#by id
function qUpdate!(network::Network, q::Real, id::Symbol)
    for element in network.elements[id]
        element.q = q
    end

    forceDensities!(network)
end

"""
Get initial lengths for form finding;
assumes units of E, A are consistent with L
"""
function initialLengths(network::Network, E::Real, A::Real)
    
    n = length(network.elements) #number of elements
    Id = I(n) #identity matrix n × n
    Em = E * Id #diagonal matrix of stiffness, E
    Am = A * Id #diagonal matrix of areas, A
    L = spdiagm(memberLengths(network)) #diagonal matrix of final lengths

    return diag((Id + (Em * Am) \ network.Q * L) \ Id)
end

"""
Initial length method for varying section properties
"""
function initialLengths(network::Network, E::Vector{<:Real}, A::Vector{<:Real})
    n = length(network.elements)

    # make sure material property vectors are the same
    @assert n == length(E) == length(A) "E and A vectors must be equal length"

    # matrix representation of components
    Id = I(n)
    Em = diagm(E)
    Am = diagm(A)
    L = spdiagm(memberLengths(network))

    return diag((Id + (Em * Am) \ network.Q * L) \ Id)
end