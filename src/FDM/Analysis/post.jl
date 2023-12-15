"""
    xyz_update!(network::Network)

Update network.xyz to reflect new equilibrium positions of nodes.
"""
function xyz_update!(network::Network)
    for index in network.N
        network.nodes[index].position = network.xyz[index, :]
    end
end

"""
    reactions!(network::Network)

Get the reaction force vectors of fixed nodes.
"""
function reactions!(network::Network)
    
    for i in network.F

        rxn = zeros(3)
        
        node = network.nodes[i]
        
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

        node.reaction = rxn
    end

end