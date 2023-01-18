abstract type AbstractNode end

"""
A node in a system
"""
mutable struct Node <: AbstractNode
    position::Vector{Float64}
    dof::Vector{Bool}
    x::Float64
    y::Float64
    z::Float64
    globalID::Vector{Int64}
    reaction::Vector{Float64}
    displacement::Vector{Float64}
    id::Union{Symbol, Nothing}

    """
    Generate a node:\\
    -position: [x,y,z] as a vector of floats\\
    -dofs: [bool, bool, bool, bool, bool, bool] as a vector of booleans where true = free DOF
    """
    function Node(position::Vector{Float64}, dofs::Vector{Bool})
        if length(position) != 3 || length(dofs) != 6
            error("Position vector must be in R³, DOFs must be length 6")
        end

        node = new(position, dofs)

        node.x, node.y, node.z = node.position

        node.id = nothing
        return node
    end

    """
    Generate a node with a common boundary condition:\\
    -position: [x,y,z] as a vector of floats\\
    -fixity: choose from :fixed, :pinned, :free, :(x/y/z)fixed, :(x/y/z)free

    """
    function Node(position::Vector{Float64}, fixity::Symbol)
        if length(position) != 3
            error("Position vector must be in R³")
        end

        dofs = copy(fixDict[fixity])

        node = new(position, dofs)
        node.x, node.y, node.z = node.position
        node.id = nothing
        return node
    end
end

mutable struct TrussNode <: AbstractNode
    position::Vector{Float64}
    dof::Vector{Bool}
    x::Float64
    y::Float64
    z::Float64
    globalID::Vector{Int64}
    reaction::Vector{Float64}
    displacement::Vector{Float64}
    id::Union{Symbol, Nothing}

    """
    Generate a node:\\
    -position: [x,y,z] as a vector of floats\\
    -dofs: [bool, bool, bool, bool, bool, bool] as a vector of booleans where true = free DOF
    """
    function TrussNode(position::Vector{Float64}, dofs::Vector{Bool})
        if length(position) != 3 || length(dofs) != 3
            error("Position vector and DOF vector must be in R³")
        end

        node = new(position, dofs)

        node.x, node.y, node.z = node.position

        node.id = nothing
        return node
    end

    """
    Generate a node with a common boundary condition:\\
    -position: [x,y,z] as a vector of floats\\
    -fixity: choose from :fixed, :pinned, :free, :(x/y/z)fixed, :(x/y/z)free

    """
    function TrussNode(position::Vector{Float64}, fixity::Symbol)
        if length(position) != 3
            error("Position vector must be in R³")
        end

        dofs = copy(fixDict[fixity][1:3])

        node = new(position, dofs)
        node.x, node.y, node.z = node.position
        node.id = nothing
        return node
    end
end

"""
Common fixity types
"""
const fixDict = Dict(:fixed => [false, false, false, false, false, false],
    :free => [true, true, true, true, true, true],
    :xfixed => [false, true, true, true, true, true],
    :yfixed => [true, false, true, true, true, true],
    :zfixed => [true, true, false, true, true, true],
    :xfree => [true, false, false, false, false, false],
    :yfree => [false, true, false, false, false, false],
    :zfree => [false, false, true, false, false, false],
    :pinned => [false, false, false, true, true, true])


const planeDict = Dict(:XY => [3, 4, 5],
    :YZ => [1, 5, 6],
    :ZX => [2, 4, 6])

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