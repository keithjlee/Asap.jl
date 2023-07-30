"""
Get initial lengths for form finding;
assumes units of E, A are consistent with L
"""
function initialLengths(network::Network, E::Real, A::Real)
    
    n = length(network.elements) #number of elements
    Id = I(n) #identity matrix n × n
    Em = E * Id #diagonal matrix of stiffness, E
    Am = A * Id #diagonal matrix of areas, A
    L = spdiagm(norm.(eachrow(network.C * network.xyz))) #diagonal matrix of final lengths

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
    L = spdiagm(norm.(eachrow(network.C * network.xyz)))

    return diag((Id + (Em * Am) \ network.Q * L) \ Id)
end

function forces(network::Network)
    return norm.(eachrow(network.C * network.xyz)) .* network.q
end