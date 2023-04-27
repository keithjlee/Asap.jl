
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
    planarize!(nodes::Vector{AbstractNode}; plane = :XY)

Restrict the dofs of a set of nodes to remain on a plane. Choose from: :XY, :YZ, :ZX
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
    fixnode!(node::AbstractNode, fixity::Symbol)

Change the DOFs of a node to a common boundary condition.
Available boundary conditions:
- :free
- :fixed
- :pinned
- :(x/y/z)free
- :(x/y/z)fixed
"""
function fixnode!(node::AbstractNode, fixity::Symbol)
    node.dof = copy(fixDict[fixity])
end

"""
Distance between two nodes
"""
function dist(n1::AbstractNode, n2::AbstractNode)
    return norm(n1.position .- n2.position)
end
