abstract type AbstractModel end

"""
Model
"""
mutable struct Model <: AbstractModel
    nodes::Vector{Node}
    elements::Vector{Element}
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
    
    function Model(nodes::Vector{Node}, elements::Vector{Element}, loads::Vector{T}) where T <: Load
        structure = new(nodes, elements, loads)
        structure.processed = false

        return structure
    end

    function Model(nodes::Vector{Node}, elements::Vector{Element})
        structure = new(nodes, elements)
        structure.loads = Vector{Load}()
        structure.processed = false

        return structure
    end
end

"""
Truss model
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
        structure.processed = false

        return structure
    end

    function TrussModel(nodes::Vector{TrussNode}, elements::Vector{TrussElement})
        structure = new(nodes, elements)
        structure.loads = Vector{NodeForce}()
        structure.processed = false

        return structure
    end
end
