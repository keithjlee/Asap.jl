"""
Reaction forces
"""
function reactions!(model::Model)
    model.reactions = zeros(model.nDOFs)
    model.reactions[model.fixedDOFs] = model.S[model.fixedDOFs, :] * model.u + model.Pf[model.fixedDOFs]
end

"""
Reaction forces
"""
function reactions!(model::TrussModel)
    model.reactions = zeros(model.nDOFs)
    model.reactions[model.fixedDOFs] = model.S[model.fixedDOFs, :] * model.u
end

"""
Populate node reactions and displacements
"""
function postprocessnodes!(model::AbstractModel)
    for node in model.nodes
        node.reaction = model.reactions[node.globalID]
        node.displacement = model.u[node.globalID]
    end
end

"""
Populate elemental forces
"""
function postprocesselements!(model::Model)
    for element in model.elements
        element.forces = element.R * (element.K * model.u[element.globalID] + element.Q)
    end
end

"""
Populate elemental forces
"""
function postprocesselements!(model::TrussModel)
    for element in model.elements
        element.forces = element.R * (element.K * model.u[element.globalID])
    end
end

"""
post-process a network
"""
function postprocess!(model::AbstractModel)
    reactions!(model)
    postprocessnodes!(model)
    postprocesselements!(model)
end