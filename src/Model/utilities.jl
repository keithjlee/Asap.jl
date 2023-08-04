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

"""
    updateDOF!(model::AbstractModel)

Update the free/fixed degrees of freedom for a model
"""
function updateDOF!(model::AbstractModel)
    model.DOFs = vcat(getproperty.(model.nodes, :dof)...)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)
end

"""
    volume(model::AbstractModel)

Get the material volume of a structural model
"""
function volume(model::AbstractModel)
    dot(getproperty.(model.elements, :length), getproperty.(getproperty.(model.elements, :section), :A))
end

function Base.copy(model::TrussModel)
    
    #new model
    nodes = Vector{TrussNode}()
    elements = Vector{TrussElement}()
    loads = Vector{NodeForce}()

    #new nodes
    for node in model.nodes
        newnode = TrussNode(copy(node.position), node.dof)
        newnode.id = node.id
        push!(nodes, newnode)
    end

    #new elements
    for element in model.elements
        newelement = TrussElement(nodes, copy(element.nodeIDs), element.section)
        newelement.id = element.id
        push!(elements, newelement)
    end

    #new loads
    for load in model.loads
        newload = NodeForce(nodes[load.node.nodeID], copy(load.value))
        newload.id = load.id
        push!(loads, newload)
    end

    model = TrussModel(nodes, elements, loads)
    solve!(model)

    return model

end