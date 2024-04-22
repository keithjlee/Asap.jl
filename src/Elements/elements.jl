abstract type AbstractElement end
abstract type FrameElement <: AbstractElement end

"""
    Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, section::Section, id::Symbol = nothing; release = :fixedfixed)
    Element(nodeStart::Node, nodeEnd::Node, section::Section, id = nothing; release = :fixedfixed)

Create a frame element with an optional `id` tag.

# Example
```julia-repl
julia> Element(nodes, [1,2], sec)
julia> Element(node1, node2, sec, :groundfloor_element)
```

# Optional argument `release` 
This property enables decoupling of nodal DOFs with respect to the end of the element.

## Available releases:
- :fixedfixed (default) - all DOFs are tied to nodes
- :fixedfree - rotational DOFs are released at end node
- :freefixed - rotational DOFs are released at start node
- :freefree - all rotational DOFs are released (truss element)
- :joist - all rotational DOFs except torsion are released

## Example

```julia-repl
julia> Element(nodes, [1,2], sec; release = :fixedfree)
```

"""
mutable struct Element <: FrameElement
    section::Section #cross section
    nodeStart::Node #start node
    nodeEnd::Node #end position
    nodeIDs::Vector{Int64} #indices of start/end nodes 
    elementID::Int64
    globalID::Vector{Int64} #element global DOFs
    loadIDs::Vector{Int64}
    length::Float64 #length of element
    release::Symbol #
    K::Matrix{Float64} # stiffness matrix in GCS
    Q::Vector{Float64} # fixed end forces in GCS
    R::Matrix{Float64} # transformation matrix
    Ψ::Float64 #roll angle
    LCS::Vector{Vector{Float64}} #local coordinate frame (X, y, z)
    forces::Vector{Float64} #elemental forces in LCS
    id::Union{Symbol, Nothing} #optional identifier

    function Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, section::Section, id = nothing; release = :fixedfixed)
        element = new(section)

        element.nodeStart, element.nodeEnd = nodes[nodeIndex]
        element.Ψ = pi/2
        element.id = id
        element.Q = zeros(12)
        element.loadIDs = Vector{Int64}()

        element.release = release

        return element
    end

    function Element(nodeStart::Node, nodeEnd::Node, section::Section, id = nothing; release = :fixedfixed)

        @assert in(release, releases) "Release not recognized; choose from: :fixedfixed, :freefixed, :fixedfree, :freefree, :joist"

        element = new(section)
        element.nodeStart = nodeStart
        element.nodeEnd = nodeEnd

        element.Ψ = pi/2
        element.id = id
        element.Q = zeros(12)
        element.loadIDs = Vector{Int64}()

        element.release = release

        return element
    end

end

"""
    BridgeElement(elementStart::Element, posStart::Float64, elementEnd::Element, posEnd::Float64, section::Section, id = nothing; release = :fixedfixed)

Create a bridge element between two frame elements. Connects from `elementStart` at a position `elementStart.length * posStart` away from `elementStart.nodeStart.position` to `elementEnd` at `elementEnd.length * posEnd` away from `elementEnd.nodeStart.position`. IE `posStart, posEnd ∈ ]0, 1[`
"""
mutable struct BridgeElement <: FrameElement
    elementStart::Element
    posStart::Float64
    elementEnd::Element
    posEnd::Float64
    section::Section
    release::Symbol
    Ψ::Float64
    elementID::Int64
    id::Union{Symbol, Nothing}
    loadIDs::Vector{Int64}

    function BridgeElement(elementStart::Element, 
            posStart::Float64, 
            elementEnd::Element, 
            posEnd::Float64, 
            section::Section,
            id = nothing;
            release = :fixedfixed)

        @assert 0 < posStart < 1 && 0 < posEnd < 1 "posStart/End must be ∈ ]0,1["
        @assert in(release, releases)

        be = new(elementStart, posStart, elementEnd, posEnd, section, release)
        be.Ψ = pi/2
        be.id = id
        be.loadIDs = Vector{Int64}()

        return be
    end
end



"""
    TrussElement(nodes::Vector{TrussNode}, nodeIndex::Vector{Int64}, section::AbstractSection, id = nothing)
    TrussElement(nodeStart::TrussNode, nodeEnd::TrussNode, section::AbstractSection, id = nothing)

Create a truss element.

# Example
```julia-repl
julia> TrussElement(nt, [1,2], sec)
TrussElement(Section(794.0, 200000.0, 77000.0, 737000.0, 737000.0, 1.47e6, 1.0), [1, 2], TrussNode([0.8879630592102802, 0.6713937498337156, 0.617463764682365], Bool[0, 0, 0], #undef, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], nothing), TrussNode([0.16046742214916832, 0.15869760269854827, 0.6247762072043447], Bool[1, 1, 1], #undef, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], nothing), #undef, 0.8900341077991537, #undef, #undef, #undef, 1.5707963267948966, #undef, nothing)
```

"""
mutable struct TrussElement <: AbstractElement
    section::AbstractSection #cross section
    nodeStart::TrussNode #start position
    nodeEnd::TrussNode #end position
    nodeIDs::Vector{Int64} #indices of start/end nodes 
    elementID::Int64
    globalID::Vector{Int64} #element global DOFs
    length::Float64 #length of element
    K::Matrix{Float64} # stiffness matrix in GCS
    R::Matrix{Float64} # transformation matrix
    forces::Vector{Float64} #elemental forces in LCS
    Ψ::Float64
    LCS::Vector{Vector{Float64}}
    id::Union{Symbol, Nothing} #optional identifier

    function TrussElement(nodes::Vector{TrussNode}, nodeIndex::Vector{Int64}, section::AbstractSection, id = nothing)
        element = new(section)

        element.nodeStart, element.nodeEnd = nodes[nodeIndex]
        element.id = id

        element.Ψ = pi/2

        return element
    end

    function TrussElement(nodeStart::TrussNode, nodeEnd::TrussNode, section::AbstractSection, id = nothing)
        element = new(section)
        element.nodeStart = nodeStart
        element.nodeEnd = nodeEnd
        element.id = id
        element.Ψ = pi/2

        return element
    end

end