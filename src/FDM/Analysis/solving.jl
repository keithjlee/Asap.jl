"""
    solve!(network::Network; reprocess = false)

Solve an FDM problem
"""
function solve!(network::Network; reprocess = false)
    if !network.processed || reprocess
        process!(network)
    end

    # solve for final free node positions
    network.xyz[network.N, :] = (network.Cn' * network.Q * network.Cn) \ (network.Pn .- network.Cn' * network.Q * network.Cf * network.xyz[network.F, :])

    # update nodal positions
    xyz_update!(network)

    #get reactions 
    reactions!(network)
end

"""
    solve(network::Network, q::Union{Vector{Int64}, Vector{Float64}})

Solve an FDM problem w/r/t a new vector of force densities
"""
function solve(network::Network, q::Union{Vector{Int64}, Vector{Float64}})

    @assert length(q) == length(network.elements) "q and elements must have equal length"

    Q = spdiagm(q)

    xyzout = copy(network.xyz)
    xyzout[network.N, :] .= (network.Cn' * Q * network.Cn) \ (network.Pn - network.Cn' * Q * network.Cf * network.xyz[network.F, :])

    return xyzout
end