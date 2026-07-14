"""
    clear_supports!(model::Model)

Remove all fixed DOFs from a model and reprocess. Will require new definitions of fixed supports before running `solve!()`
"""
function clear_supports!(model::Model)

    for node in model.nodes
        fixnode!(node, :free)
        node.id = :free
    end

    process!(model)

end

"""
    element_connectivity(model::Model)

Get the Cel = [n_element × n_element] adjacency matrix for a given model, where Cel[i, j] = 1 if element i shares a node with element j, and is 0 otherwise.
"""
function element_connectivity(model::Model)

    #element-node connectivity [ne × nn]
    C = connectivity(model)

    #initialize element-element connectivity matrix [ne × ne]
    n_elements = length(model.elements)
    Cel = spzeros(Int64, n_elements, n_elements)


    for irow in axes(C, 1)
        i_nodes = findall(C[irow, :] .!= 0)

        for icol in i_nodes
            col = C[:, icol]

            i_connected = findall(col .!= 0)

            for i_element in i_connected
                if i_element == irow
                    continue
                else
                    Cel[irow, i_element] = 1
                end
            end
        end
    end

    return Cel

end
