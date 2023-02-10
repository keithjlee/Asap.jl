
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

"""
Set the DOF of a node to a common type
"""
function fixnode!(node::Node, fixity::Symbol)
    node.dof = copy(fixDict[fixity])
end

"""
Distance between two nodes
"""
function dist(n1::AbstractNode, n2::AbstractNode)
    return norm(n1.position .- n2.position)
end
