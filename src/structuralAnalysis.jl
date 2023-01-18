"""
populate global DOF indices of nodes and elements
"""
function populateDOF!(model::Model)

    model.DOFs = vcat([node.dof for node in model.nodes]...)
    model.nDOFs = length(model.DOFs)
    model.nNodes = length(model.nodes)
    model.nElements = length(model.elements)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)

    n_dof = 6
    for (i, node) in enumerate(model.nodes)
        node.globalID = i * n_dof - (n_dof - 1) .+ collect(0:n_dof-  1)
    end

    for element in model.elements
        idStart = model.nodes[element.nodeIDs[1]].globalID
        idEnd = model.nodes[element.nodeIDs[2]].globalID

        element.globalID = [idStart; idEnd]
    end
end

"""
populate global DOF indices of nodes and elements
"""
function populateDOF!(model::TrussModel)

    model.DOFs = vcat([node.dof for node in model.nodes]...)
    model.nDOFs = length(model.DOFs)
    model.nNodes = length(model.nodes)
    model.nElements = length(model.elements)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)

    n_dof = 3
    for (i, node) in enumerate(model.nodes)
        node.globalID = i * n_dof - (n_dof - 1) .+ collect(0:n_dof-  1)
    end

    for element in model.elements
        idStart = model.nodes[element.nodeIDs[1]].globalID
        idEnd = model.nodes[element.nodeIDs[2]].globalID

        element.globalID = [idStart; idEnd]
    end
end

"""
Process elements: get transformation matrix and global elemental stiffness matrix
"""
function processElements!(model::Model)
    for element in model.elements
        element.Q = zeros(12) # reset Qf
        element.R = R(element)
        element.LCS = lcs(element, element.Ψ)
        makeK!(element)
    end
end

"""
Process elements: get transformation matrix and global elemental stiffness matrix
"""
function processElements!(model::TrussModel)
    for element in model.elements
        element.R = R(element)
        element.LCS = lcs(element, element.Ψ)
        makeK!(element)
    end
end

"""
Populate a node force load
"""
function populateLoad!(model::Model, load::NodeForce)
    idx = load.node.globalID[1:3]
    model.P[idx] += load.value
end

"""
Populate a node force load
"""
function populateLoad!(model::TrussModel, load::NodeForce)
    idx = load.node.globalID
    model.P[idx] += load.value
end


function populateLoad!(P::Vector{Float64}, load::NodeForce)
    idx = load.node.globalID[1:3]
    P[idx] += load.value
end

"""
Populate a node moment load
"""
function populateLoad!(model::Model, load::NodeMoment)
    idx = load.node.globalID[4:6]
    model.P[idx] += load.value
end

function populateLoad!(P::Vector{Float64}, load::NodeMoment)
    idx = load.node.globalID[4:6]
    P[idx] += load.value
end

"""
Populate element loads
"""
function populateLoad!(model::Model, load::ElementLoad)
    idx = load.element.globalID
    load.element.Q += load.element.R' * Q(load)
    model.Pf[idx] += load.element.Q
end

function populateLoad!(P::Vector{Float64}, load::ElementLoad)
    idx = load.element.globalID
    P[idx] += load.element.R' * Q(load)
end

"""
Create load vectors P and Pf
"""
function populateLoads!(model::Model)
    #initialize
    model.P = zeros(model.nDOFs)
    model.Pf = zeros(model.nDOFs)

    #create load vectors
    for load in model.loads
        populateLoad!(model, load)
    end
end

"""
Create load vectors P and Pf
"""
function populateLoads!(model::TrussModel)
    #initialize
    model.P = zeros(model.nDOFs)

    #create load vectors
    for load in model.loads
        populateLoad!(model, load)
    end
end

"""
create load vector F = P-Pf
"""
function createF(model::Model, loads::Vector{Load})

    P = zeros(model.nDOFs)
    Pf = zeros(model.nDOFs)

    for load in loads
        if typeof(load) <: ElementLoad
            populateLoad!(Pf, load)
        else
            populateLoad!(P, load)
        end
    end

    return P - Pf
end

"""
create load vector F = P-Pf
"""
function createF(model::TrussModel, loads::Vector{NodeForce})

    P = zeros(model.nDOFs)

    for load in loads
        populateLoad!(P, load)
    end

    return P
end



"""
Create global stiffness matrix
"""
function globalS!(model::Union{Model, TrussModel})
    S = spzeros(model.nDOFs, model.nDOFs)

    for element in model.elements
        idx = element.globalID
        S[idx, idx] .+= element.K
    end

    model.S = Symmetric(S)
end

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
function postprocessnodes!(model::Union{Model, TrussModel})
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
process a network
"""
function process!(model::Union{Model, TrussModel})

    #global DOF 
    populateDOF!(model)

    #elements 
    processElements!(model)

    #loads
    populateLoads!(model)

    #stiffness matrix
    globalS!(model)

    #processing finished
    model.processed = true
end

"""
post-process a network
"""
function postprocess!(model::Union{Model, TrussModel})
    reactions!(model)
    postprocessnodes!(model)
    postprocesselements!(model)
end

"""
Solve a network
"""
function solve!(model::Model; reprocess = false)

    if !model.processed || reprocess
        process!(model)
    end

    # reduce scope of problem and solve
    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    U = model.S[idx, idx] \ F

    #compliance
    model.compliance = U' * F

    #full DOF displacement vector
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    # post process
    postprocess!(model)
end

"""
Solve a network
"""
function solve!(model::TrussModel; reprocess = false)

    if !model.processed || reprocess
        process!(model)
    end

    # reduce scope of problem and solve
    idx = model.freeDOFs
    U = model.S[idx, idx] \ model.P[idx]

    #compliance
    model.compliance = U' * model.P[idx]

    #full DOF displacement vector
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    # post process
    postprocess!(model)
end

"""
Solve a network and return displacement vector
"""
function solve(model::Model)
    idx = model.freeDOFs

    F = model.P[idx] - model.Pf[idx]

    U = model.S[idx, idx] \ F

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
Solve a network and return displacement vector
"""
function solve(model::TrussModel)
    idx = model.freeDOFs

    F = model.P[idx]

    U = model.S[idx, idx] \ F

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
Solve a network with a new load vector
"""
function solve(model::Model, F::Vector{Float64})
    if length(F) != model.nDOFs
        error("Length of F must equal total number of DOFs")
    end

    idx = model.freeDOFs

    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
Solve a network with a new load vector
"""
function solve(model::TrussModel, F::Vector{Float64})

    if length(F) != model.nDOFs
        error("Length of F must equal total number of DOFs")
    end

    idx = model.freeDOFs

    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
Solve a network with a new set of loads
"""
function solve(model::Model, L::Vector{Load})
    
    F = createF(model, L)
    
    idx = model.freeDOFs

    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
Solve a network with a new set of loads
"""
function solve(model::TrussModel, L::Vector{NodeForce})
    
    F = createF(model, L)
    
    idx = model.freeDOFs

    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end