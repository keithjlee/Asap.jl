"""
applied to network
"""
function xyz_update!(network::Network)
    for index in network.N
        network.nodes[index].position = network.xyz[index, :]
    end
end

function reactions!(network::Network)
    
    for i in network.F

        rxn = zeros(3)
        
        node = network.nodes[i]
        
        iels = findall(network.C[:, i] .!= 0)

        for index in iels
            e = network.elements[index]
            eforce = vector(e) * force(e)

            if network.C[index, i] < 0
                rxn .+= eforce
            else
                rxn .-= eforce
            end
        end

        node.reaction = rxn
    end

end

"""
∑|F|L
"""
function FL(network::Network)
    return sum(norm.(eachrow(network.C * network.xyz)).^2 .* network.q)
end