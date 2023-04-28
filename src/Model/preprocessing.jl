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

    # assign an id to node, extract global DOF index
    for (i, node) in enumerate(model.nodes)
        node.nodeID = i
        node.globalID = i * n_dof - (n_dof - 1) .+ collect(0:n_dof-  1)
    end

    #assign an id to load, store load id into relevant node/element
    for (i, load) in enumerate(model.loads)
        load.loadID = i

        assign!(load)
    end

    # assign an id to element, get associated node IDs, extract global DOF
    for (i, element) in enumerate(model.elements)
        element.elementID = i
        element.nodeIDs = [element.nodeStart.nodeID, element.nodeEnd.nodeID]

        idStart = element.nodeStart.globalID
        idEnd = element.nodeEnd.globalID

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
        node.nodeID = i
        node.globalID = i * n_dof - (n_dof - 1) .+ collect(0:n_dof - 1)
    end

    for (i, load) in enumerate(model.loads)
        load.loadID = i

        assign!(load)
    end

    for (i, element) in enumerate(model.elements)

        element.elementID = i
        element.nodeIDs = [element.nodeStart.nodeID, element.nodeEnd.nodeID]

        idStart = element.nodeStart.globalID
        idEnd = element.nodeEnd.globalID

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
        element.length = dist(element.nodeStart, element.nodeEnd)
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
        element.length = dist(element.nodeStart, element.nodeEnd)
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
function globalS!(model::Model)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()

    for element in model.elements

        idx = element.globalID

        @inbounds for i = 1:12
            @inbounds for j = 1:12
                k = element.K[i,j]
                push!(I, idx[i])
                push!(J, idx[j])
                push!(V, k)
                
            end
        end

    end

    model.S = sparse(I, J, V)
end

"""
Create global stiffness matrix
"""
function globalS!(model::TrussModel)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()

    for element in model.elements

        idx = element.globalID

        @inbounds for i = 1:6
            @inbounds for j = 1:6
                k = element.K[i,j]
                push!(I, idx[i])
                push!(J, idx[j])
                push!(V, k)
                
            end
        end

    end

    model.S = sparse(I, J, V)
end

"""
process a network
"""
function process!(model::AbstractModel)

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

