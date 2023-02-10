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