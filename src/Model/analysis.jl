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
    solve(model::AbstractModel)

Solve a model and directly return the global displacement vector.
"""
function solve(model::Model)
    idx = model.freeDOFs

    F = model.P[idx] - model.Pf[idx]

    U = model.S[idx, idx] \ F

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

function solve(model::TrussModel)
    idx = model.freeDOFs

    F = model.P[idx]

    U = model.S[idx, idx] \ F

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
    solve(model::AbstractModel, F::Vector{Float64})

Find the displacement of the structural model w/r/t a new external force vector. 
"""
function solve(model::Model, F::Vector{Float64})

    @assert length(F) == model.nDOFs "Length of F must equal total number of DOFs"

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

    @assert length(F) == model.nDOFs "Length of F must equal total number of DOFs"

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