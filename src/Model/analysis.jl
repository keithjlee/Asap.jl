"""
    process!(model::Model)

Process a structural model: add linkages between nodes and elements, determine DOF orders, generate the load vectors P, Pf, and assemble the global stiffness matrix, S.
"""
function process!(model::Model)

    make_ids!(model)

    if any(typeof.(model.elements) .== BridgeElement)
        processBridge!(model)
    else
        process_elements!(model)
    end

    #global DOF 
    populate_DOF_indices!(model)

    #loads
    populate_loads!(model)

    #stiffness matrix
    create_S!(model)

    #processing finished
    model.processed = true
end

"""
    process!(model::TrussModel)

Process a structural truss model: add linkages between nodes and elements, determine DOF orders, generate the load vector P, and assemble the global stiffness matrix, S.
"""
function process!(model::TrussModel)

    make_ids!(model)

    process_elements!(model)

    #global DOF 
    populate_DOF_indices!(model)

    #loads
    populate_loads!(model)

    #stiffness matrix
    create_S!(model)

    #processing finished
    model.processed = true
end

"""
    solve!(model::Model; reprocess = false)

Solve for the nodal displacements of a structural model. `reprocess = true` reevaluates all node/element properties and reassembles the global stiffness matrix.
"""
function solve!(model::Model; reprocess = false)

    if !model.processed || reprocess
        # clear existing load associations
        for node in model.nodes
            empty!(node.loadIDs)
        end

        for element in model.elements
            empty!(element.loadIDs)
            element.Q = zero(element.Q)
        end
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
    post_process!(model)
end

"""
    solve!(model::TrussModel; reprocess = false)

Solve for the nodal displacements of a structural truss model. `reprocess = true` reevaluates all node/element properties and reassembles the global stiffness matrix.
"""
function solve!(model::TrussModel; reprocess = false)

    if !model.processed || reprocess
        # clear existing load associations
        for node in model.nodes
            empty!(node.loadIDs)
        end
    
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
    post_process!(model)
end

"""
    solve(model::Model, L::Vector{AbstractLoad})

Return the displacement vector under a given load set L.
"""
function solve(model::Model, L::Vector{<:AbstractLoad})

   model.processed  || process!(model)
    
    F = create_F(model, L)
    
    idx = model.freeDOFs

    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
    solve!(model::Model, L::Vector{AbstractLoad})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::Model, L::Vector{<:AbstractLoad})

    model.loads = L

    for node in model.nodes
        empty!(node.loadIDs)
    end

    for element in model.elements
        empty!(element.loadIDs)
    end

    process!(model)
    solve!(model)

    # post process
    post_process!(model)
end

"""
    solve(model::TrussModel, L::Vector{NodeForce})

Return the displacement vector to a new set of loads.
"""
function solve(model::TrussModel, L::Vector{NodeForce})

    !model.processed || process!(model)
    
    F = create_F(model, L)
    
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

    model.loads = L

    for node in model.nodes
        empty!(node.loadIDs)
    end

    process!(model)
    solve!(model)

    # post process
    post_process!(model)
end
