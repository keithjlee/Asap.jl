abstract type AbstractElement end

"""
An element
"""
mutable struct Element <: AbstractElement
    section::Section #cross section
    nodeIDs::Vector{Int64} #indices of start/end nodes 
    posStart::Vector{Float64} #start position
    posEnd::Vector{Float64} #end position
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

    """
    Base constructor of a frame element:\\
    -nodes: vector of nodes used in model\\
    -nodeIndex: [startIndex, endIndex] of element w/r/t nodes\\
    -section: element Section
    """
    function Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, section::Section)
        element = new(section, nodeIndex)

        element.posStart = nodes[nodeIndex[1]].position
        element.posEnd = nodes[nodeIndex[2]].position
        element.length = norm(element.posEnd .- element.posStart)
        element.Ψ = pi/2
        # element.LCS = lcs(element, element.Ψ)
        element.id = nothing
        element.Q = zeros(12)

        element.release = :fixedfixed

        return element
    end

    """
    Element constructor with prescribed end releases:\\
    -nodes: vector of nodes used in model\\
    -nodeIndex: [startIndex, endIndex] of element w/r/t nodes\\
    -section: element Section
    -release: choose from :fixedfixed, :freefixed, :fixedfree, :freefree

    """
    function Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, section::Section, release::Symbol)

        if !in(release, releases)
            error("Release not recognized; choose from: :fixedfixed, :freefixed, :fixedfree, :freefree, :truss")
        end

        element = new(section, nodeIndex)

        element.posStart = nodes[nodeIndex[1]].position
        element.posEnd = nodes[nodeIndex[2]].position
        element.length = norm(element.posEnd .- element.posStart)
        element.Ψ = pi/2
        # element.LCS = lcs(element, element.Ψ)
        element.id = nothing
        element.Q = zeros(12)

        element.release = release

        return element
    end

end


"""
A truss element
"""
mutable struct TrussElement <: AbstractElement
    section::AbstractSection #cross section
    nodeIDs::Vector{Int64} #indices of start/end nodes 
    posStart::Vector{Float64} #start position
    posEnd::Vector{Float64} #end position
    globalID::Vector{Int64} #element global DOFs
    length::Float64 #length of element
    K::Matrix{Float64} # stiffness matrix in GCS
    R::Matrix{Float64} # transformation matrix
    forces::Vector{Float64} #elemental forces in LCS
    Ψ::Float64
    LCS::Vector{Vector{Float64}}
    id::Union{Symbol, Nothing} #optional identifier

    """
    Base constructor of a frame element:\\
    -nodes: vector of nodes used in model\\
    -nodeIndex: [startIndex, endIndex] of element w/r/t nodes\\
    -section: element Section
    """
    function TrussElement(nodes::Vector{TrussNode}, nodeIndex::Vector{Int64}, section::AbstractSection)
        element = new(section, nodeIndex)

        element.posStart = nodes[nodeIndex[1]].position
        element.posEnd = nodes[nodeIndex[2]].position
        element.length = norm(element.posEnd .- element.posStart)
        element.id = nothing

        element.Ψ = pi/2

        return element
    end

end