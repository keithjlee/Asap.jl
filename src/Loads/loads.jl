"""
Define loads applied to structure/elements
"""
abstract type AbstractLoad end
abstract type NodeLoad <: AbstractLoad end
abstract type ElementLoad <: AbstractLoad end

"""
    NodeForce(node::AbstractNode, value::Vector{Float64})
    NodeForce(nodes::Vector{<:AbstractNode}, index::Integer, value::Vector{Float64})

A force vector [Fx, Fy, Fz] in the global coordinate system applied to a node.
"""
mutable struct NodeForce <: NodeLoad
    node::Union{Node, TrussNode}
    value::Vector{Float64}
    loadID::Int64
    id::Symbol
    
    function NodeForce(node::AbstractNode, value::Vector{Float64}, id::Symbol = :force)

        @assert length(value) == 3 "load vector must be in R³ (GCS)"

        force = new(node, value, 0, id)

        return force
    end

    function NodeForce(nodes::Vector{<:AbstractNode}, index::Integer, value::Vector{Float64}, id::Symbol = :force)

        NodeForce(nodes[index], value, id)

        return force
    end
end

"""
    NodeMoment(node::Node, value::Vector{Float64})

A moment vector [Mx, My, Mz] in the global coordinate system applied to a node with rotational DOFs.
"""
mutable struct NodeMoment <: NodeLoad
    node::Node
    value::Vector{Float64}
    loadID::Int64
    id::Symbol
    
    function NodeMoment(node::Node, value::Vector{Float64}, id::Symbol = :moment)

        @assert length(value) == 3 "Moment vector must be in R³ (GCS)"

        force = new(node, value, 0, id)

        return force
    end
end

"""
    LineLoad(element::Element, value::Vector{Float64})

A distributed line load [wx, wy, wz] in (force/length) applied along an element in the global coordinate system.
"""
mutable struct LineLoad <: ElementLoad
    element::FrameElement
    value::Vector{Float64}
    loadID::Int64
    id::Symbol

    function LineLoad(element::T, value::Vector{Float64}, id::Symbol = :lineload) where T <: FrameElement

        @assert length(value) == 3 "load vector must be in R³ (GCS)"

        force = new(element, value, 0, id)

        return force
    end
end

"""
    GravityLoad(element::Element, factor::Float64)

A gravity load (negative global Z) applied along a member. 

Generates distributed load w = element.section.A * element.section.ρ * factor, where factor should be the appropriate acceleration due to gravity.
"""
mutable struct GravityLoad <: ElementLoad
    element::FrameElement
    factor::Float64
    loadID::Int64
    id::Symbol

    function GravityLoad(element::FrameElement, factor::Float64, id::Symbol = :gravityload)
        force = new(element, factor, 0, id)
        return force
    end
end


"""
    PointLoad(element::Element, position::Float64, value::Vector{Float64})

A point load [Px, Py, Pz] applied in the global coordinate system at a distance `position` × `element.length` from the starting node.
"""
mutable struct PointLoad <: ElementLoad
    element::FrameElement
    position::Float64
    value::Vector{Float64}
    loadID::Int64
    id::Symbol

    function PointLoad(element::FrameElement, position::Float64, value::Vector{Float64}, id::Symbol = :pointload)
        @assert 0 < position < 1 "position must be ∈ ]0, 1["
        @assert length(value) == 3 "load vector must be in R³ (GCS)"

        force = new(element, position, value, 0, id)
        return force
    end
end