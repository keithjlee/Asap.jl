"""
    solve!(network::Network; reprocess = false)

Solve an FDM network for its equilibrium geometry and update the node
positions in place (free axes only — fixed coordinates are boundary data).

FDM equilibrium is SEPARABLE per coordinate: the force-density matrix
`D = CᵀQC` depends only on topology and force densities, never on
positions. Uniform networks (every node fully free or fully fixed) solve
all three coordinates in one multi-RHS system; networks with per-axis
fixity (e.g. plan prescribed, height free) solve each axis against its own
free/fixed partition.
"""
function solve!(network::Network; reprocess = false)
    if !network.processed || reprocess
        process!(network)
    end

    if network.mixed
        for a in 1:3
            Na = network.Naxis[a]
            isempty(Na) && continue
            Fa = network.Faxis[a]
            Cn = network.C[:, Na]
            Cf = network.C[:, Fa]
            network.xyz[Na, a] = (Cn' * network.Q * Cn) \
                (network.P[Na, a] .- Cn' * network.Q * Cf * network.xyz[Fa, a])
        end
    else
        # uniform: one factorization, three right-hand sides
        network.xyz[network.N, :] = (network.Cn' * network.Q * network.Cn) \
            (network.Pn .- network.Cn' * network.Q * network.Cf * network.xyz[network.F, :])
    end

    # update nodal positions
    update_xyz!(network)

    #get reactions
    reactions!(network)
end

"""
    solve(network::Network, q::Vector{<:Real}) -> xyz

Equilibrium positions for a NEW vector of force densities without mutating
the network (same per-axis semantics as [`solve!`](@ref)).
"""
function solve(network::Network, q::Vector{<:Real})
    @assert length(q) == length(network.elements) "q and elements must have equal length"

    Q = spdiagm(q)
    xyzout = copy(network.xyz)

    if network.mixed
        for a in 1:3
            Na = network.Naxis[a]
            isempty(Na) && continue
            Fa = network.Faxis[a]
            Cn = network.C[:, Na]
            Cf = network.C[:, Fa]
            xyzout[Na, a] = (Cn' * Q * Cn) \
                (network.P[Na, a] .- Cn' * Q * Cf * network.xyz[Fa, a])
        end
    else
        xyzout[network.N, :] .= (network.Cn' * Q * network.Cn) \
            (network.Pn - network.Cn' * Q * network.Cf * network.xyz[network.F, :])
    end

    return xyzout
end
