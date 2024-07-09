abstract type AbstractElement end
abstract type FrameElement{R} <: AbstractElement end

# Element release type parameters
abstract type Release end
struct FixedFixed <: Release end
struct FixedFree <: Release end
struct FreeFixed <: Release end
struct FreeFree <: Release end
struct Joist <: Release end

const _ReleaseDict = Dict(
    :fixedfixed => FixedFixed,
    :fixedfree => FixedFree,
    :freefixed => FreeFixed,
    :freefree => FreeFree,
    :joist => Joist
)

"""
    Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, section::Section, id::Symbol = nothing; release = :fixedfixed)
    Element(nodeStart::Node, nodeEnd::Node, section::Section, id = :element; release = :fixedfixed)

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
mutable struct Element{R<:Release} <: FrameElement{R}
    section::Section #cross section
    nodeStart::Node #start node
    nodeEnd::Node #end position
    elementID::Int64
    globalID::Vector{Int64} #element global DOFs
    length::Float64 #length of element
    K::Matrix{Float64} # stiffness matrix in GCS
    Q::Vector{Float64} # fixed end forces in GCS
    R::Matrix{Float64} # transformation matrix
    Ψ::Float64 #roll angle
    LCS::Vector{Vector{Float64}} #local coordinate frame (X, y, z)
    forces::Vector{Float64} #elemental forces in LCS
    id::Symbol #optional identifier

    function Element(nodeStart::Node, nodeEnd::Node, section::Section, id = :element; release = :fixedfixed)

        @assert in(release, keys(_ReleaseDict)) "Release not recognized; choose from: :fixedfixed, :freefixed, :fixedfree, :freefree, :joist"

        element = new{_ReleaseDict[release]}(
            section,
            nodeStart,
            nodeEnd,
            0,
            Vector{Int64}(undef, 12),
            0.0,
            zeros(12,12),
            zeros(12),
            zeros(12,12),
            pi/2,
            repeat([zeros(3)], 3),
            zeros(12),
            id
        )

        return element
    end

end

"""
    BridgeElement(elementStart::Element, posStart::Float64, elementEnd::Element, posEnd::Float64, section::Section, id = :element; release = :fixedfixed)

Create a bridge element between two frame elements. Connects from `elementStart` at a position `elementStart.length * posStart` away from `elementStart.nodeStart.position` to `elementEnd` at `elementEnd.length * posEnd` away from `elementEnd.nodeStart.position`. IE `posStart, posEnd ∈ ]0, 1[`
"""
mutable struct BridgeElement{R<:Release} <: FrameElement{R}
    elementStart::Element
    posStart::Float64
    elementEnd::Element
    posEnd::Float64
    section::Section
    release::Symbol
    Ψ::Float64
    elementID::Int64
    id::Union{Symbol, Nothing}

    function BridgeElement(elementStart::Element, 
            posStart::Float64, 
            elementEnd::Element, 
            posEnd::Float64, 
            section::Section,
            id = :element;
            release = :fixedfixed)

        @assert 0 < posStart < 1 && 0 < posEnd < 1 "posStart/End must be ∈ ]0,1["
        @assert in(release, keys(_ReleaseDict))

        be = new{_ReleaseDict[release]}(
            elementStart, 
            posStart, 
            elementEnd, 
            posEnd, 
            section, 
            release, 
            pi/2, 
            0, 
            id
        )

        return be
    end
end



"""
    TrussElement(nodes::Vector{TrussNode}, nodeIndex::Vector{Int64}, section::AbstractSection, id = :element)
    TrussElement(nodeStart::TrussNode, nodeEnd::TrussNode, section::AbstractSection, id = :element)

Create a truss element.

# Example
```julia-repl
julia> TrussElement(nt, [1,2], sec)
TrussElement(Section(794.0, 200000.0, 77000.0, 737000.0, 737000.0, 1.47e6, 1.0), [1, 2], TrussNode([0.8879630592102802, 0.6713937498337156, 0.617463764682365], Bool[0, 0, 0], #undef, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], nothing), TrussNode([0.16046742214916832, 0.15869760269854827, 0.6247762072043447], Bool[1, 1, 1], #undef, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], nothing), #undef, 0.8900341077991537, #undef, #undef, #undef, 1.5707963267948966, #undef, nothing)
```

"""
mutable struct TrussElement <: AbstractElement
    section::Union{TrussSection,Section} #cross section
    nodeStart::TrussNode #start position
    nodeEnd::TrussNode #end position
    elementID::Int64
    globalID::Vector{Int64} #element global DOFs
    length::Float64 #length of element
    K::Matrix{Float64} # stiffness matrix in GCS
    R::Matrix{Float64} # transformation matrix
    forces::Vector{Float64} #elemental forces in LCS
    Ψ::Float64
    LCS::Vector{Vector{Float64}}
    id::Union{Symbol, Nothing} #optional identifier

    function TrussElement(nodeStart::TrussNode, nodeEnd::TrussNode, section::AbstractSection, id = :element)
        
        element = new(
            section,
            nodeStart,
            nodeEnd,
            0,
            Vector{Int64}(undef, 6),
            0.0,
            zeros(2,2),
            zeros(2,6),
            zeros(6),
            pi/2,
            repeat([zeros(3)], 3),
            id
        )

        return element
    end

end