"""
    NodalSpring{T<:Real}

An elastic (spring) support acting at a node — **applicative data**: the
spring references its node the way a load references its target, and lives
in the model's spring list. Nodes never know about their springs; multiple
springs on one node compose additively.

Each of the node's six global DOFs gets an independent spring stiffness:
translational stiffnesses [force/length] on (Tx, Ty, Tz), rotational
stiffnesses [force·length/rad] on (Rx, Ry, Rz). Zero entries mean no spring
on that DOF. During assembly the stiffnesses are added to the global
stiffness diagonal; a sprung DOF is marked *active* (it participates in the
solve even if nothing else touches it), and its spring reaction is
recovered as `−k·u` in post-processing.

Use a `NodalSpring` for elastic foundations, soil springs, flexible
supports, and partial restraints — anywhere a rigid `fixity` support is too
stiff an idealization.

# Fields
- `node::Node{T}`: the supported node
- `stiffness::SVector{6,T}`: spring stiffness per global DOF,
  ordered (Tx, Ty, Tz, Rx, Ry, Rz); entries ≥ 0
- `id::Symbol`: user tag

# Constructors
    NodalSpring(node, stiffness::AbstractVector, id = :spring)   # 6 entries
    NodalSpring(node, k_translation::Real, id = :spring)         # uniform Tx=Ty=Tz=k

# Examples
```julia-repl
julia> soil = NodalSpring(base, [0.0, 0.0, 5e4, 0.0, 0.0, 0.0], :soil)  # vertical only

julia> pad = NodalSpring(base, 1e5)   # equal translational springs, no rotational
```
"""
struct NodalSpring{T<:Real}
    node::Node{T}
    stiffness::SVector{6,T}
    id::Symbol

    function NodalSpring(node::Node{T}, stiffness::AbstractVector{<:Real},
        id::Symbol=:spring) where {T}
        @assert length(stiffness) == 6 "spring stiffness must have 6 entries (Tx, Ty, Tz, Rx, Ry, Rz)"
        @assert all(k -> k >= 0, stiffness) "spring stiffnesses must be ≥ 0"
        return new{T}(node, SVector{6,T}(stiffness), id)
    end
end

NodalSpring(node::Node{T}, k::Real, id::Symbol=:spring) where {T} =
    NodalSpring(node, SVector{6,T}(k, k, k, 0, 0, 0), id)
