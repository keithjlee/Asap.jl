"""
    update_xyz!(network::Network)

Copy the solved coordinate matrix back onto the node positions. Fixed
coordinates are boundary data and unchanged; every free axis takes its
equilibrium value.
"""
function update_xyz!(network::Network)
    for (i, node) in enumerate(network.nodes)
        any(node.fixity) || continue
        for a in 1:3
            node.fixity[a] && (node.position[a] = network.xyz[i, a])
        end
    end
end

"""
    reactions!(network::Network)

Anchor forces on every FIXED node axis: the sum of connected member forces
along that axis. Components on free axes are zero — equilibrium already
holds there, so there is nothing for the support to supply.
"""
function reactions!(network::Network)
    for (i, node) in enumerate(network.nodes)
        all(node.fixity) && continue        # fully free: no reaction

        rxn = zeros(3)
        iels = findall(network.C[:, i] .!= 0)
        for index in iels
            e = network.elements[index]
            eforce = local_x(e) * force(e)
            if network.C[index, i] < 0
                rxn .+= eforce
            else
                rxn .-= eforce
            end
        end

        for a in 1:3
            node.fixity[a] && (rxn[a] = 0.0)   # free axis: not a reaction
        end
        node.reaction = rxn
    end
end
