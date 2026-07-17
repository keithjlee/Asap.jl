"""
    initial_lengths(network::Network, E::Real, A::Real)
    initial_lengths(network::Network, E::Vector, A::Vector)

Required unstressed member lengths so that, at the solved geometry, members
of axial stiffness `E·A` carry exactly the form-found forces `q·L`
[length]. Solve the network first.
"""
function initial_lengths(network::Network, E::Real, A::Real)
    n = Base.length(network.elements)
    return initial_lengths(network, fill(E, n), fill(A, n))
end

function initial_lengths(network::Network, E::Vector{<:Real}, A::Vector{<:Real})
    n = Base.length(network.elements)
    @assert n == Base.length(E) == Base.length(A) "E and A vectors must be equal length"

    Id = I(n)
    Em = diagm(E)
    Am = diagm(A)
    q = [el.q for el in network.elements]
    L = spdiagm(Base.length.(network.elements))    # lengths at current geometry

    return diag((Id + (Em * Am) \ (spdiagm(q) * L)) \ Id)
end

"""
    forces(network::Network) -> Vector

Member axial forces `q·L` at the current geometry [force].
"""
forces(network::Network) = [force(el) for el in network.elements]

"""
    update_q!(network::Network, q)
    update_q!(network::Network, q::Vector, indices::Vector{Int})

Set force densities and re-solve — a convenience for form-finding sweeps
(assembly reads `q` fresh, so this is just assignment plus `solve!`).
"""
function update_q!(network::Network{T}, q::Real) where {T}
    for el in network.elements
        el.q = T(q)
    end
    solve!(network)
end

function update_q!(network::Network{T}, q::Vector{<:Real}) where {T}
    @assert Base.length(q) == Base.length(network.elements)
    for (el, qi) in zip(network.elements, q)
        el.q = T(qi)
    end
    solve!(network)
end

function update_q!(network::Network{T}, q::Vector{<:Real}, indices::Vector{Int}) where {T}
    @assert Base.length(q) == Base.length(indices)
    for (i, qi) in zip(indices, q)
        network.elements[i].q = T(qi)
    end
    solve!(network)
end
