"""
process a network
"""
function process!(model::Model)

    if any(typeof.(model.elements) .== BridgeElement)
        processBridge!(model)
    else
        processElements!(model)
    end

    #global DOF 
    populateDOF!(model)

    #loads
    populateLoads!(model)

    #stiffness matrix
    globalS!(model)

    #processing finished
    model.processed = true
end

function process!(model::TrussModel)

    processElements!(model)

    #global DOF 
    populateDOF!(model)

    #loads
    populateLoads!(model)

    #stiffness matrix
    globalS!(model)

    #processing finished
    model.processed = true
end

"""
    solve!(model::AbstractModel; reprocess = false)

Perform a structural analysis. `reprocess = true` re-generates the stiffness matrix and resets all saved solutions.
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
    solve(model::Model, L::Vector{Load})

Return the displacement vector to a new set of loads
"""
function solve(model::Model, L::Vector{<:Load})

    @assert model.processed "Model must be already processed"
    
    F = createF(model, L)
    
    idx = model.freeDOFs

    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
    solve!(model::Model, L::Vector{Load})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::Model, L::Vector{<:Load})

    @assert model.processed "Model must be already processed"

    model.loads = L

    # clear existing load associations
    for node in model.nodes
        empty!(node.loadIDs)
    end

    for element in model.elements
        empty!(element.loadIDs)
    end

    # assign new load associations
    for (i, load) in enumerate(model.loads)
        load.loadID = i
        assign!(load)
    end

    # create P, Pf matrix
    populateLoads!(model)
    
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
    solve(model::TrussModel, L::Vector{NodeForce})

Return the displacement vector to a new set of loads
"""
function solve(model::TrussModel, L::Vector{NodeForce})

    !model.processed && process!(model)
    
    F = createF(model, L)
    
    idx = model.freeDOFs

    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
    solve!(model::TrussModel, L::Vector{NodeForce})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::TrussModel, L::Vector{NodeForce})

    @assert model.processed "Model must be already processed"

    model.loads = L

    # clear existing load associations
    for node in model.nodes
        empty!(node.loadIDs)
    end

    # assign new load associations
    for (i, load) in enumerate(model.loads)
        load.loadID = i
        assign!(load)
    end

    # create P, Pf matrix
    populateLoads!(model)
    
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
