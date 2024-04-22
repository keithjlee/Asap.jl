
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

import Base: -
-(node1::AbstractNode, node2::AbstractNode) = node1.position - node2.position

"""
    planarize!(nodes, plane = :XY)

Restrict the DOFs of a set of nodes to a given global plane. E.g., when analyzing a 2D structure.

# Arguments
- `nodes::Vector{<:AbstractNode}` vector of nodes to restrict
- `plane::Symbol = :XY` plane to restrict DOFs. Defaults to the XY plane, can be:
    - :XY
    - :XZ
    - :YZ
"""
function planarize!(nodes::Vector{Node}, plane = :XY)
    idx = planeDict[plane]

    for node in nodes
        node.dof[idx] .= false
    end
end

function planarize!(nodes::Vector{TrussNode}, plane = :XY)
    idx = planeDict[plane][1]

    for node in nodes
        node.dof[idx] = false
    end
end

"""
    fixnode!(node::AbstractNode, fixity::Symbol)

Fix the DOFs of a node to a common boundary condition.

# Arguments
- `node::AbstractNode`node to modify
- `fixity::Symbol` boundary condition to apply. Available boundary conditions:
    - :free
    - :fixed
    - :pinned
    - :(x/y/z)free
    - :(x/y/z)fixed
"""
function fixnode!(node::Node, fixity::Symbol)
    node.dof = copy(fixDict[fixity])
end

function fixnode!(node::TrussNode, fixity::Symbol)
    node.dof = copy((fixDict[fixity])[1:3])
end