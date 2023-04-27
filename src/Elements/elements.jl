abstract type AbstractElement end

"""
    Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, section::Section)

Instantiate a frame element with a given section that connects two nodes.

# Example
```julia-repl
julia> Element(nodes, [1,2], sec)
Element(Section(794.0, 200000.0, 77000.0, 737000.0, 737000.0, 1.47e6, 1.0), [1, 2], Node([0.0, 0.0, 0.0], Bool[0, 0, 0, 0, 0, 0], #undef, #undef, #undef, nothing), Node([5000.0, 500.0, 5200.0], Bool[0, 0, 0, 1, 1, 1], #undef, #undef, #undef, nothing), #undef, 7231.18247591637, :fixedfixed, #undef, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], #undef, 1.5707963267948966, #undef, #undef, nothing)
```

-----------------------------

    Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, section::Section, release::Symbol)

Instantiate a frame element with a given section that connects two nodes with a given end release.
Available releases:
- :fixedfixed (default)
- :fixedfree
- :freefixed
- :freefree

# Example

```julia-repl
julia> Element(nodes, [1,2], sec, :fixedfree)
Element(Section(794.0, 200000.0, 77000.0, 737000.0, 737000.0, 1.47e6, 1.0), [1, 2], Node([0.0, 0.0, 0.0], Bool[0, 0, 0, 0, 0, 0], #undef, #undef, #undef, nothing), Node([5000.0, 500.0, 5200.0], Bool[0, 0, 0, 1, 1, 1], #undef, #undef, #undef, nothing), #undef, 7231.18247591637, :fixedfree, #undef, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], #undef, 1.5707963267948966, #undef, #undef, nothing)
```

"""
mutable struct Element <: AbstractElement
    section::Section #cross section
    nodeIDs::Vector{Int64} #indices of start/end nodes 
    nodeStart::Node #start node
    nodeEnd::Node #end position
    globalID::Vector{Int64} #element global DOFs
    length::Float64 #length of element
    release::Symbol #
    K::Matrix{Float64} # stiffness matrix in GCS
    Q::Vector{Float64} # fixed end forces in GCS
    R::Matrix{Float64} # transformation matrix
    Ψ::Float64 #roll angle
    LCS::Vector{Vector{Float64}} #local coordinate frame (X, y, z)
    forces::Vector{Float64} #elemental forces in LCS
    id::Union{Symbol, Nothing} #optional identifier

    function Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, section::Section)
        element = new(section, nodeIndex)

        element.nodeStart, element.nodeEnd = nodes[nodeIndex]
        element.length = dist(element.nodeStart, element.nodeEnd)
        element.Ψ = pi/2
        element.id = nothing
        element.Q = zeros(12)

        element.release = :fixedfixed

        return element
    end

    function Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, section::Section, release::Symbol)

        if !in(release, releases)
            error("Release not recognized; choose from: :fixedfixed, :freefixed, :fixedfree, :freefree, :truss")
        end

        element = new(section, nodeIndex)

        element.nodeStart, element.nodeEnd = nodes[nodeIndex]
        element.length = dist(element.nodeStart, element.nodeEnd)
        element.Ψ = pi/2
        # element.LCS = lcs(element, element.Ψ)
        element.id = nothing
        element.Q = zeros(12)

        element.release = release

        return element
    end

end


"""
    TrussElement(nodes::Vector{TrussNode}, nodeIndex::Vector{Int64}, section::AbstractSection)

Instantiate a truss element with a given section that connects two truss nodes.

# Example
```julia-repl
julia> TrussElement(nt, [1,2], sec)
TrussElement(Section(794.0, 200000.0, 77000.0, 737000.0, 737000.0, 1.47e6, 1.0), [1, 2], TrussNode([0.8879630592102802, 0.6713937498337156, 0.617463764682365], Bool[0, 0, 0], #undef, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], nothing), TrussNode([0.16046742214916832, 0.15869760269854827, 0.6247762072043447], Bool[1, 1, 1], #undef, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], nothing), #undef, 0.8900341077991537, #undef, #undef, #undef, 1.5707963267948966, #undef, nothing)
```

"""
mutable struct TrussElement <: AbstractElement
    section::AbstractSection #cross section
    nodeIDs::Vector{Int64} #indices of start/end nodes 
    nodeStart::TrussNode #start position
    nodeEnd::TrussNode #end position
    globalID::Vector{Int64} #element global DOFs
    length::Float64 #length of element
    K::Matrix{Float64} # stiffness matrix in GCS
    R::Matrix{Float64} # transformation matrix
    forces::Vector{Float64} #elemental forces in LCS
    Ψ::Float64
    LCS::Vector{Vector{Float64}}
    id::Union{Symbol, Nothing} #optional identifier

    function TrussElement(nodes::Vector{TrussNode}, nodeIndex::Vector{Int64}, section::AbstractSection)
        element = new(section, nodeIndex)

        element.nodeStart, element.nodeEnd = nodes[nodeIndex]
        element.length = dist(element.nodeStart, element.nodeEnd)
        element.id = nothing

        element.Ψ = pi/2

        return element
    end

end