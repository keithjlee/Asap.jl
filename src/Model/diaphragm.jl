"""
    RigidDiaphragm{T}

A rigid in-plane diaphragm: the in-plane translations and the rotation
about the plane normal of every SLAVE node follow a MASTER node as a rigid
body — the conventional floor model for building lateral analysis (a floor
slab is far stiffer in its plane than the walls and frames it ties
together). Out-of-plane translation and the two out-of-plane rotations of
each slave stay independent, so floor bending, column axial forces and
wall overturning are unaffected.

For a diaphragm normal to global axis `a` with in-plane axes `(b, c)`
(Z: (X, Y); X: (Y, Z); Y: (Z, X)) and offsets `r = x_slave − x_master`:

    u_b,s = u_b,m − θ_a,m · r_c
    u_c,s = u_c,m + θ_a,m · r_b
    θ_a,s = θ_a,m

Imposed exactly by transformation (`u = C·û`), not by penalty springs — see
`docs/shells/DERIVATION.md` §5.

# Fields
- `nodes::Vector{Node{T}}`: slave nodes (the floor's nodes)
- `master::Node{T}`: the node whose 3 in-plane DOFs carry the diaphragm;
  apply story forces and torsional moments here
- `axis::Int`: global axis of the plane normal (1 = X, 2 = Y, 3 = Z)
- `id::Symbol`

# Constructor
    RigidDiaphragm(nodes; master = nothing, axis = :Z, id = :diaphragm)

With `master = nothing` a master node is created at the slaves' centroid:
its 3 in-plane DOFs free, its 3 out-of-plane DOFs fixed (it is a kinematic
node, not a structural one). The model adds it to its node list.

# Example
```julia-repl
julia> floor2 = RigidDiaphragm(level2_nodes; id = :L2);

julia> model = Model(nodes, elements, [NodeForce(floor2.master, [Fx, 0, 0])];
                     diaphragms = [floor2])
```

Slave nodes must not be supported in the diaphragm's in-plane DOFs, and a
node may belong to only one diaphragm. Gradients through the pure path
treat the master–slave offsets as constant (frozen at the last `solve!`).
"""
struct RigidDiaphragm{T}
    nodes::Vector{Node{T}}
    master::Node{T}
    axis::Int
    id::Symbol

    function RigidDiaphragm(nodes::AbstractVector{Node{T}}; master::Union{Nothing,Node{T}} = nothing,
        axis::Symbol = :Z, id::Symbol = :diaphragm) where {T}
        a = axis === :X ? 1 : axis === :Y ? 2 : axis === :Z ? 3 :
            throw(ArgumentError("diaphragm axis must be :X, :Y or :Z"))
        length(nodes) >= 2 || throw(ArgumentError("a diaphragm needs at least two slave nodes"))
        allunique(nodes) || throw(ArgumentError("diaphragm slave nodes must be distinct"))
        if master === nothing
            c = sum(n -> n.position, nodes) / length(nodes)
            fix = trues(6)
            fix[a] = false                          # out-of-plane translation
            fix[3 + mod1(a + 1, 3)] = false         # the two out-of-plane rotations
            fix[3 + mod1(a + 2, 3)] = false
            master = Node(collect(c), collect(fix), id)
        end
        any(n -> n === master, nodes) && throw(ArgumentError("the master cannot also be a slave"))
        return new{T}(collect(nodes), master, a, id)
    end
end

# in-plane axes (b, c) of a diaphragm normal to axis a
@inline _inplane_axes(a::Int) = (mod1(a + 1, 3), mod1(a + 2, 3))
