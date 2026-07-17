"""
    FDMnode(x, y, z, dof::Bool, id = :node)
    FDMnode(pos, dof::Bool, id = :node)
    FDMnode(pos, fixity::AbstractVector{Bool}, id = :node)

A node in an FDM network: a spatial position and a PER-AXIS degree of
freedom. `fixity` follows the same convention as the core `Node`
(`true` = free, `false` = fixed), so a node may be free in some axes and
fixed in others — e.g. `[false, false, true]` prescribes the plan position
and lets the height find equilibrium (vault form-finding). The `Bool`
constructors set all three axes at once (`true` = fully free).

# Fields
- `position::Vector{Float64}`: the [x, y, z] position [length] — updated in
  place by `solve!` for every free axis (the network geometry IS the result)
- `fixity::SVector{3,Bool}`: per-axis freedom (`true` = free)
- `id::Symbol`: identifier (optional)
- `nodeID::Int`: global index (internal, assigned by `process!`)
- `reaction::Vector{Float64}`: anchor force [force] on the FIXED axes after
  `solve!` (free-axis components are zero — equilibrium holds there)
"""
mutable struct FDMnode
    position::Vector{Float64}
    fixity::SVector{3,Bool}
    id::Symbol
    nodeID::Int64
    reaction::Vector{Float64}

    #empty constructor
    function FDMnode()
        return new()
    end

    function FDMnode(x::Real, y::Real, z::Real, dof::Bool, id = :node)
        return new(Float64.([x, y, z]), SVector(dof, dof, dof), id)
    end

    function FDMnode(pos::Vector{<:Real}, dof::Bool, id = :node)
        @assert length(pos) == 3 "pos should be length 3"
        return new(Float64.(pos), SVector(dof, dof, dof), id)
    end

    function FDMnode(pos::Vector{<:Real}, fixity::AbstractVector{Bool}, id = :node)
        @assert length(pos) == 3 "pos should be length 3"
        @assert length(fixity) == 3 "fixity should be length 3 (true = free)"
        return new(Float64.(pos), SVector{3,Bool}(fixity), id)
    end
end
include("utilities.jl")
