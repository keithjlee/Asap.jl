abstract type AbstractModel end

function makeids!(nodes::Vector{<:AbstractNode})
    @inbounds for (i, node) in enumerate(nodes)
        node.nodeID = i
    end
end

function makeids!(elements::Vector{<:AbstractElement})
    @inbounds for (i, element) in enumerate(elements)
        element.elementID = i
    end
end

function makeids!(loads::Vector{<:Load})
    @inbounds for (i, load) in enumerate(loads)
        load.loadID = i
    end
end

"""
    Model(nodes::Vector{Node}, elements::Vector{Element}, loads::Vector{Load})

Create a complete structural model ready for analysis.
"""
mutable struct Model <: AbstractModel
    nodes::Vector{Node}
    elements::Vector{FrameElement}
    loads::Vector{Load}
    nNodes::Int64
    nElements::Int64
    DOFs::Vector{Bool} #vector of DOFs
    nDOFs::Int64
    freeDOFs::Vector{Int64} #free DOF indices
    fixedDOFs::Vector{Int64}
    S::SparseMatrixCSC{Float64,Int64} # global stiffness
    P::Vector{Float64} # external loads
    Pf::Vector{Float64} # element end forces
    u::Vector{Float64} # nodal displacements
    reactions::Vector{Float64} # reaction forces
    compliance::Float64 #structural compliance
    processed::Bool
    
    function Model(nodes::Vector{Node}, elements::Vector{<:FrameElement}, loads::Vector{<:Load})
        structure = new(nodes, elements, loads)

        makeids!(structure.nodes)
        makeids!(structure.elements)
        makeids!(structure.loads)

        structure.processed = false

        return structure
    end

    # function Model(nodes::Vector{Node}, elements::Vector{<:FrameElement})
    #     structure = new(nodes, elements)
    #     structure.loads = Vector{Load}()

    #     makeids!(structure.nodes)
    #     makeids!(structure.elements)

    #     structure.processed = false

    #     return structure
    # end
end

"""
    TrussModel(nodes::Vector{TrussNode}, elements::Vector{TrussElement}, loads::Vector{NodeForce})

Create a complete structural model ready for analysis.

"""
mutable struct TrussModel <: AbstractModel
    nodes::Vector{TrussNode}
    elements::Vector{TrussElement}
    loads::Vector{NodeForce}
    nNodes::Int64
    nElements::Int64
    DOFs::Vector{Bool} #vector of DOFs
    nDOFs::Int64
    freeDOFs::Vector{Int64} #free DOF indices
    fixedDOFs::Vector{Int64}
    S::SparseMatrixCSC{Float64,Int64} # global stiffness
    P::Vector{Float64} # external loads
    u::Vector{Float64} # nodal displacements
    reactions::Vector{Float64} # reaction forces
    compliance::Float64 #structural compliance
    processed::Bool
    
    function TrussModel(nodes::Vector{TrussNode}, elements::Vector{TrussElement}, loads::Vector{NodeForce})
        structure = new(nodes, elements, loads)

        makeids!(structure.nodes)
        makeids!(structure.elements)
        makeids!(structure.loads)

        structure.processed = false

        return structure
    end

    # function TrussModel(nodes::Vector{TrussNode}, elements::Vector{TrussElement})
    #     structure = new(nodes, elements)
    #     structure.loads = Vector{NodeForce}()

    #     makeids!(structure.nodes)
    #     makeids!(structure.elements)

    #     structure.processed = false

    #     return structure
    # end
end
