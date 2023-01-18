# Custom indexing
function Base.getindex(elements::Vector{Element}, i::Symbol)
    return [element for element in elements if element.id == i]
end

function Base.findall(elements::Vector{Element}, i::Symbol)
    return findall([element.id == i for element in elements])
end

function Base.getindex(elements::Vector{TrussElement}, i::Symbol)
    return [element for element in elements if element.id == i]
end

function Base.findall(elements::Vector{TrussElement}, i::Symbol)
    return findall([element.id == i for element in elements])
end

function Base.getindex(nodes::Vector{Node}, i::Symbol)
    return [node for node in nodes if node.id == i]
end

function Base.findall(nodes::Vector{Node}, i::Symbol)
    return findall([node.id == i for node in nodes])
end

function Base.getindex(nodes::Vector{TrussNode}, i::Symbol)
    return [node for node in nodes if node.id == i]
end

function Base.findall(nodes::Vector{TrussNode}, i::Symbol)
    return findall([node.id == i for node in nodes])
end

function Base.getindex(loads::Vector{Load}, i::Symbol)
    return [load for load in loads if load.id == i]
end

function Base.findall(loads::Vector{Load}, i::Symbol)
    return findall([load.id == i for load in loads])
end

"""
Fix nodes to a plane. All DOF (including rotations) that would cause a non-planar deformation are set to false.\\
-nodes: a vector of nodes to planarize
-plane: choose from :XY, :YZ, :ZX
"""
function planarize!(nodes::Vector{Node}; plane = :XY)
    idx = planeDict[plane]

    for node in nodes
        node.dof[idx] .= false
    end
end

function planarize!(nodes::Vector{TrussNode}; plane = :XY)
    idx = planeDict[plane][1]

    for node in nodes
        node.dof[idx] = false
    end
end

function planarize!(model::AbstractModel; plane = :XY)
    planarize!(model.nodes; plane = plane)
    if plane == :XY
        for element in model.elements
            element.Ψ = 0.
        end
    end
end