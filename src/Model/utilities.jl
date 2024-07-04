"""
    planarize!(model::AbstractModel, plane = :XY)

Fix all nodal DOFs to remain on plane = `plane`
"""
function planarize!(model::AbstractModel, plane = :XY)
    planarize!(model.nodes, plane)
    if plane == :XY
        for element in model.elements
            element.Ψ = 0.
        end
    end
end

"""
    connectivity(model::AbstractModel)

Get the [nₑ × nₙ] sparse matrix where C[i, j] = -1 if element i starts at node j, and C[i,j] = 1 if element i ends at node j, and 0 otherwise.
"""
function connectivity(model::AbstractModel)
    I = vcat([[i, i] for i = 1:model.nElements]...)
    J = vcat([nodeids(e) for e in model.elements]...)
    V = repeat([-1, 1], model.nElements)

    return sparse(I, J, V)
end

"""
    node_positions(model::AbstractModel)

Generate the [nₙ × 3] node position matrix
"""
function node_positions(model::AbstractModel)
    return vcat([node.position' for node in model.nodes]...)
end

export nodePositions
nodePositions(model::AbstractModel) = node_positions(model)

"""
    update_DOF!(model::AbstractModel)

Update the free/fixed degrees of freedom for a model
"""
function update_DOF!(model::AbstractModel)
    model.DOFs = vcat(getproperty.(model.nodes, :dof)...)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)
end

export updateDOF!
updateDOF!(model::AbstractModel) = update_DOF!(model)

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
        newnode = TrussNode(deepcopy(node.position), node.dof)
        newnode.id = node.id
        push!(nodes, newnode)
    end

    #new elements
    for element in model.elements
        newelement = TrussElement(nodes, deepcopy(element.nodeIDs), element.section)
        newelement.id = element.id
        push!(elements, newelement)
    end

    #new loads
    for load in model.loads
        newload = NodeForce(nodes[load.node.nodeID], deepcopy(load.value))
        newload.id = load.id
        push!(loads, newload)
    end

    model = TrussModel(nodes, elements, loads)
    solve!(model)

    return model

end