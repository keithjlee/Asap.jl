"""
    populate_DOF_indices!(model::Model)

Populate all indices, references, and other information in a model.
"""
function populate_DOF_indices!(model::Model)

    model.DOFs = vcat(getproperty.(model.nodes, :dof)...)
    model.nDOFs = length(model.DOFs)
    model.nNodes = length(model.nodes)
    model.nElements = length(model.elements)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)

    n_dof = 6
    dofset = collect(0:n_dof - 1)

    # assign an id to node, extract global DOF index
    for (i, node) in enumerate(model.nodes)
        node.globalID = i * n_dof - (n_dof - 1) .+ dofset
    end

    # assign an id to element, get associated node IDs, extract global DOF
    for element in model.elements
        element.globalID = [element.nodeStart.globalID; element.nodeEnd.globalID]
    end

end

"""
    populate_DOF_indices!(model::TrussModel)

Populate all indices, references, and other information in a truss model.
"""
function populate_DOF_indices!(model::TrussModel)

    model.DOFs = vcat(getproperty.(model.nodes, :dof)...)
    model.nDOFs = length(model.DOFs)
    model.nNodes = length(model.nodes)
    model.nElements = length(model.elements)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)

    n_dof = 3
    for (i, node) in enumerate(model.nodes)
        node.globalID = i * n_dof - (n_dof - 1) .+ collect(0:n_dof - 1)
    end

    for element in model.elements
        element.globalID = [element.nodeStart.globalID; element.nodeEnd.globalID]
    end

end

"""
    process_elements!(model::Model)
    
Populate the transformation matrix and global elemental stiffness matrix of the elements in a model.
"""
function process_elements!(model::Model)
    for element in model.elements
        fill!(element.Q, 0.0)
        lcs!(element, element.Ψ)
        R!(element)
        length!(element)
        global_K!(element)
    end
end

"""
    process_elements!(elements::Vector{<:FrameElement})
    
Populate the transformation matrix and global elemental stiffness matrix of the elements in a vector of elements.
"""
function process_elements!(elements::Vector{T}) where {T<:FrameElement}
    for element in elements
        fill!(element.Q, 0.0)
        lcs!(element, element.Ψ)
        R!(element)
        length!(element)
        global_K!(element)
    end
end

"""
    process_elements!(model::TrussModel)
    
Populate the transformation matrix and global elemental stiffness matrix of the elements in a truss model.
"""
function process_elements!(model::TrussModel)
    for element in model.elements
        lcs!(element, element.Ψ)
        R!(element)
        length!(element)
        global_K!(element)
    end
end

"""
    populate_load!(model::AbstractModel, load::NodeForce)

Populate the global load vector `model.P` with a nodal force.
"""
function populate_load!(model::AbstractModel, load::NodeForce)
    idx = load.node.globalID[1:3]
    model.P[idx] += load.value
end

"""
    populate_load!(P::Vector{Float64}, load::NodeForce)

Populate the global load vector P` with a nodal force.
"""
function populate_load!(P::Vector{Float64}, load::NodeForce)
    idx = load.node.globalID[1:3]
    P[idx] += load.value
end

"""
    populate_load!(model::Model, load::NodeMoment)

Populate the global load vector `model.P` with a nodal moment.
"""
function populate_load!(model::Model, load::NodeMoment)
    idx = load.node.globalID[4:6]
    model.P[idx] += load.value
end


"""
    populate_load!(model::Model, load::NodeMoment)

Populate the global load vector `model.P` with a nodal moment.
"""
function populate_load!(P::Vector{Float64}, load::NodeMoment)
    idx = load.node.globalID[4:6]
    P[idx] += load.value
end

"""
    populate_load!(model::Model, load::ElementLoad)

Generate the fixed-end force vector `Q` for a given load, and populate the global fixed-end force vector `Pf`.
"""
function populate_load!(model::Model, load::ElementLoad)
    idx = load.element.globalID
    fixed_end_forces_GCS = load.element.R' * q(load)
    load.element.Q += fixed_end_forces_GCS
    model.Pf[idx] += fixed_end_forces_GCS
end

"""
    populate_load!(P::Vector{Float64}, load::ElementLoad)

Populate an external fixed-end force vector Pf with respect the an elemental load `load`
"""
function populate_load!(Pf::Vector{Float64}, load::ElementLoad)
    idx = load.element.globalID
    Pf[idx] += load.element.R' * q(load)
end

"""
    populate_loads!(model::Model)

Generate the nodal force vectors `model.P` (external) and `model.Pf` (fixed-end)
"""
function populate_loads!(model::Model)

    #clear load vectors
    fill!(model.P, 0.0)
    fill!(model.Pf, 0.0)

    #create load vectors
    for load in model.loads
        populate_load!(model, load)
    end
end

"""
    populate_loads!(model::TrussModel)

Generate the nodal force vectors `model.P`
"""
function populate_loads!(model::TrussModel)
    
    #clear load vectors
    fill!(model.P, 0.0)

    #create load vectors
    for load in model.loads
        populate_load!(model, load)
    end
end

"""
    create_F(model::TrussModel, loads::Vector{NodeForce})

create load vector F = P - Pf
"""
function create_F(model::Model, loads::Vector{T}) where {T<:AbstractLoad}

    P = zeros(model.nDOFs)
    Pf = zeros(model.nDOFs)

    for load in loads
        if typeof(load) <: ElementLoad
            populate_load!(Pf, load)
        else
            populate_load!(P, load)
        end
    end

    return P - Pf
end

"""
    create_F(model::TrussModel, loads::Vector{NodeForce})

create load vector F = P
"""
function create_F(model::TrussModel, loads::Vector{NodeForce})

    P = zeros(model.nDOFs)

    for load in loads
        populate_load!(P, load)
    end

    return P
end

"""
    create_S!(model::Model)

Assemble the global stiffness matrix S.
"""
function create_S!(model::Model)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()

    for element in model.elements

        idx = element.globalID

        for i = 1:12
            for j = 1:12
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
    create_S!(model::TrussModel)

Assemble the global stiffness matrix S.
"""
function create_S!(model::TrussModel)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()

    for element in model.elements

        idx = element.globalID

        for i = 1:6
            for j = 1:6
                k = element.K[i,j]
                push!(I, idx[i])
                push!(J, idx[j])
                push!(V, k)
                
            end
        end

    end

    model.S = sparse(I, J, V)
end