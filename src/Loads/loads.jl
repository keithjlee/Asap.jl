"""
Define loads applied to structure/elements
"""
abstract type Load end
abstract type NodeLoad <: Load end
abstract type ElementLoad <: Load end

"""
A force applied to a node in global XYZ axes
"""
mutable struct NodeForce <: NodeLoad
    node::Union{Node, TrussNode}
    value::Vector{Float64}
    id::Union{Symbol, Nothing}
    
    function NodeForce(node::Union{Node, TrussNode}, value::Vector{Float64})

        if length(value) != 3
            error("Load value must be a vector in R³")
        end

        force = new(node, value)
        force.id = nothing

        return force
    end
end

"""
A moment applied to a node in global XYZ axes
"""
mutable struct NodeMoment <: NodeLoad
    node::Node
    value::Vector{Float64}
    id::Union{Symbol, Nothing}
    
    function NodeMoment(node::Node, value::Vector{Float64})

        if length(value) != 3
            error("Load value must be a vector in R³")
        end

        force = new(node, value)
        force.id = nothing

        return force
    end
end

"""
A distributed line load along an element with respect to global XYZ axes.\\
Generates w = norm(value) [force/distance] in the direction of value.
"""
mutable struct LineLoad <: ElementLoad
    element::Element
    value::Vector{Float64}
    id::Union{Symbol, Nothing}

    function LineLoad(element::Element, value::Vector{Float64})

        if length(value) != 3
            error("Load value must be a vector in R³")
        end

        force = new(element, value)
        force.id = nothing
        return force
    end
end

"""
A gravity load (negative global Z) applied along a member.\\
Generates distributed load w = element.section.A * element.section.ρ * factor\\
Where factor should be the appropriate acceleration due to gravity
"""
mutable struct GravityLoad <: ElementLoad
    element::Element
    factor::Float64
    id::Union{Symbol, Nothing}

    function GravityLoad(element::Element, factor::Float64)
        force = new(element, factor)
        force.id = nothing
        return force
    end
end


"""
A point load applied at a fractional point along element.\\
Generates a load vector P = [Px, Py, Pz] at a distance L = element.length × position from the starting node
"""
mutable struct PointLoad <: ElementLoad
    element::Element
    position::Float64
    value::Vector{Float64}
    id::Union{Symbol, Nothing}

    function PointLoad(element::Element, position::Float64, value::Vector{Float64})
        if !(0.0 < position < 1.0)
            error("position must be > 0 and < 1")
        end

        if length(value) != 3
            error("Load value must be a vector in R³")
        end

        force = new(element, position, value)
        force.id = nothing
        return force
    end
end