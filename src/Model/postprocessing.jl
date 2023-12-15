"""
    reactions!(model::AbstractModel)

Populate external reaction forces in `model.reactions`
"""
function reactions!(model::Model)
    model.reactions = zeros(model.nDOFs)
    model.reactions[model.fixedDOFs] = model.S[model.fixedDOFs, :] * model.u + model.Pf[model.fixedDOFs]
end

function reactions!(model::TrussModel)
    model.reactions = zeros(model.nDOFs)
    model.reactions[model.fixedDOFs] = model.S[model.fixedDOFs, :] * model.u
end

"""
    post_process_nodes!(model::AbstractModel)

Populate nodal reaction and displacement fields
"""
function post_process_nodes!(model::AbstractModel)
    for node in model.nodes
        node.reaction = model.reactions[node.globalID]
        node.displacement = model.u[node.globalID]
    end
end

"""
    post_process_elements!(model::AbstractModel)

Populate elemental LCS force vectors in `element.forces`
"""
function post_process_elements!(model::Model)
    for element in model.elements
        element.forces = element.R * (element.K * model.u[element.globalID] + element.Q)
    end
end

function post_process_elements!(model::TrussModel)
    for element in model.elements
        element.forces = element.R * (element.K * model.u[element.globalID])
    end
end

"""
    post_process!(model::AbstractModel)

Post process a model after solving for displacements using `solve!(model)`
"""
function post_process!(model::AbstractModel)
    reactions!(model)
    post_process_nodes!(model)
    post_process_elements!(model)
end