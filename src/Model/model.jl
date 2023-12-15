abstract type AbstractModel end

function make_ids!(nodes::Vector{<:AbstractNode})
    @inbounds for (i, node) in enumerate(nodes)
        node.nodeID = i
    end
end

function make_ids!(elements::Vector{TrussElement})
    @inbounds for (i, element) in enumerate(elements)
        element.elementID = i
    end
end

function make_ids!(elements::Vector{<:FrameElement})
    @inbounds for (i, element) in enumerate(elements)
        element.elementID = i
    end
end

function make_ids!(loads::Vector{<:AbstractLoad})
    @inbounds for (i, load) in enumerate(loads)
        load.loadID = i
    end
end

function make_ids!(model::AbstractModel)

    make_ids!(model.nodes)
    make_ids!(model.elements)
    make_ids!(model.loads)

end

"""
    Model(nodes::Vector{Node}, elements::Vector{Element}, loads::Vector{AbstractLoad})

Create a complete structural model ready for analysis.
"""
mutable struct Model <: AbstractModel
    nodes::Vector{Node}
    elements::Vector{FrameElement}
    loads::Vector{AbstractLoad}
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
    tol::Float64
    processed::Bool
    
    function Model(nodes::Vector{Node}, elements::Vector{<:FrameElement}, loads::Vector{<:AbstractLoad})
        structure = new(nodes, elements, loads)

        structure.processed = false
        structure.tol = 1e-3

        return structure
    end

    # function Model(nodes::Vector{Node}, elements::Vector{<:FrameElement})
    #     structure = new(nodes, elements)
    #     structure.loads = Vector{AbstractLoad}()

    #     make_ids!(structure.nodes)
    #     make_ids!(structure.elements)

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
    tol::Float64
    processed::Bool
    
    function TrussModel(nodes::Vector{TrussNode}, elements::Vector{TrussElement}, loads::Vector{NodeForce})
        structure = new(nodes, elements, loads)

        structure.processed = false
        structure.tol = 1e-3

        return structure
    end

    # function TrussModel(nodes::Vector{TrussNode}, elements::Vector{TrussElement})
    #     structure = new(nodes, elements)
    #     structure.loads = Vector{NodeForce}()

    #     make_ids!(structure.nodes)
    #     make_ids!(structure.elements)

    #     structure.processed = false

    #     return structure
    # end
end
