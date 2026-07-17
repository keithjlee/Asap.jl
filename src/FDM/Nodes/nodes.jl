"""
    FDMnode{T}

A node in an FDM network: a spatial position and a PER-AXIS degree of
freedom. `fixity` follows the same convention as the core `Node`
(`true` = free, `false` = fixed), so a node may be free in some axes and
fixed in others — e.g. `[false, false, true]` prescribes the plan position
and lets the height find equilibrium (vault form-finding).

Form-finding note: `solve!` updates `position` in place on every free axis —
the network geometry IS the analysis result. Reactions are NOT stored here;
query them with `reaction(network.results, node)`.

# Fields
- `position::SVector{3,T}`: the [x, y, z] position [length]
- `fixity::SVector{3,Bool}`: per-axis freedom (`true` = free)
- `id::Symbol`: identifier (optional)
- `index::Int`: global index (internal, assigned by `process!`)

# Constructors
    FDMnode(x, y, z, dof::Bool, id = :node)
    FDMnode(pos, dof::Bool, id = :node)               # all three axes at once
    FDMnode(pos, fixity::AbstractVector{Bool}, id = :node)
"""
mutable struct FDMnode{T<:Real}
    position::SVector{3,T}
    fixity::SVector{3,Bool}
    id::Symbol
    index::Int

    function FDMnode(pos::AbstractVector{<:Real}, fixity::AbstractVector{Bool}, id::Symbol=:node)
        @assert Base.length(pos) == 3 "pos should be length 3"
        @assert Base.length(fixity) == 3 "fixity should be length 3 (true = free)"
        T = float(promote_type(eltype(pos), Float64))
        return new{T}(SVector{3,T}(pos), SVector{3,Bool}(fixity), id, 0)
    end
end

FDMnode(pos::AbstractVector{<:Real}, dof::Bool, id::Symbol=:node) =
    FDMnode(pos, SVector(dof, dof, dof), id)
FDMnode(x::Real, y::Real, z::Real, dof::Bool, id::Symbol=:node) =
    FDMnode([x, y, z], dof, id)

include("utilities.jl")
