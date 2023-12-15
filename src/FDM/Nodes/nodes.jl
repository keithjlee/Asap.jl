"""
    FDMnode(x::Real, y::Real, z::Real, dof::Bool, id = nothing)
    FDMnode(pos::Vector{<:Real}, dof::Bool, id = nothing)

A node in an FDM network defined by its spatial position [x, y, z] and degree of freedom (free = `true`, fixed = `false`)

# Fields
- position::Vector{Float64}: the [x,y,z] position of the node
- dof::Bool: true = free; false = fixed
- id::Symbol: identifier (optional)
- nodeID::Integer: global indentifier (internal)
- reaction::Vector{Float64}: the [x,y,z] reaction forces if node is fixed
"""
mutable struct FDMnode
    position::Vector{Float64}
    dof::Bool # true = free; false = fixed
    id::Union{Symbol, Nothing}
    nodeID::Integer
    reaction::Vector{Float64}

    #empty constructor
    function FDMnode()
        return new()
    end

    # individual coordinate basis
    function FDMnode(x::Real, y::Real, z::Real, dof::Bool, id = nothing)
        return new(Float64.([x, y, z]), dof, id)
    end

    # # using a vector to represent position
    # function FDMnode(pos::Vector{Float64}, dof::Bool, id = nothing)
    #     @assert length(pos) == 3 "pos should be length 3"
        
    #     return FDMnode(pos, dof, id)
    # end

    function FDMnode(pos::Vector{<:Real}, dof::Bool, id = nothing)
        @assert length(pos) == 3 "pos should be length 3"
        
        return new(Float64.(pos), dof, id)
    end

end
include("utilities.jl")