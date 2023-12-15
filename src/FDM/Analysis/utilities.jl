"""
    initial_lengths(network::Network, E::Real, A::Real)

Get the required unstressed length of elements given a material stiffness E and cross sectional area A.
"""
function initial_lengths(network::Network, E::Real, A::Real)
    
    n = length(network.elements) #number of elements
    Id = I(n) #identity matrix n × n
    Em = E * Id #diagonal matrix of stiffness, E
    Am = A * Id #diagonal matrix of areas, A
    L = spdiagm(norm.(eachrow(network.C * network.xyz))) #diagonal matrix of final lengths

    return diag((Id + (Em * Am) \ network.Q * L) \ Id)
end

"""
    initial_lengths(network::Network, E::Vector{<:Real}, A::Vector{<:Real})

Get the required unstressed length of elements given an element-wise material stiffness vector E and cross sectional area vector A.
"""
function initial_lengths(network::Network, E::Vector{<:Real}, A::Vector{<:Real})
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

"""
    forces(network::Network)

Return the axial forces in a network.
"""
function forces(network::Network)
    return norm.(eachrow(network.C * network.xyz)) .* network.q
end

"""
    update_q!(network::Network, q::Float64)

Set all force density values to `q` and re-solve.
"""
function update_q!(network::Network, q::Float64)
    setfield!.(network.elements, :q, q)
    force_densities!(network)

    solve!(network)
end

"""
    update_q!(network::Network, q::Vector{Float64})

Update all force densities with `q` and re-solve.
"""
function update_q!(network::Network, q::Vector{Float64})
    @assert length(q) == length(network.elements)

    setfield!.(network.elements, :q, q)
    force_densities!(network)

    solve!(network)
end

"""
    update_q!(network::Network, q::Vector{Float64}, indices::Vector{Int64}))

Update all force densities of elements at `indices` with values in `q` and resolve.
"""
function update_q!(network::Network, q::Vector{Float64}, indices::Vector{Int64})
    @assert length(q) == length(indices)

    setfield!.(network.elements[indices], :q, q)
    force_densities!(network)

    solve!(network)
end