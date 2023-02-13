"""
Fix DOFs of all nodes to become planar
"""
function planarize!(model::AbstractModel; plane = :XY)
    planarize!(model.nodes; plane = plane)
    if plane == :XY
        for element in model.elements
            element.Ψ = 0.
        end
    end
end

"""
Generate the element-node connectivity matrix
"""
function connectivity(model::AbstractModel)
    I = vcat([[i, i] for i = 1:model.nElements]...)
    J = vcat([e.nodeIDs for e in model.elements]...)
    V = repeat([-1, 1], model.nElements)

    return sparse(I, J, V)
end

"""
Generate the (nₙ × 3) node position matrix
"""
function nodePositions(model::AbstractModel)
    return vcat([node.position' for node in model.nodes]...)
end