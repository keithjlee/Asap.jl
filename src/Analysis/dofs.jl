"""
    DofPartition

The classification of every global degree of freedom into exactly one of
three states — the backbone of the analysis layer:

- **free**: *active* (some element or spring couples stiffness to it) and
  not fixed by a support → an unknown of the solve
- **fixed**: fixed by a nodal support (`fixity = false`) → a reaction slot
- **inactive**: free of supports but touched by no stiffness → excluded
  from the system entirely. This is what removes truss-node rotations and
  fully released end rotations *structurally*: no zero-stiffness rows, no
  singular modes, no regularization.

Global DOF numbering: node `i` owns slots `6(i−1)+1 … 6i` ordered
(Tx, Ty, Tz, Rx, Ry, Rz); internal (non-nodal) element DOFs follow after
all nodal slots.

# Fields
- `n_global::Int`: total DOF slots (nodal + internal)
- `free::Vector{Int}`: global indices of solve unknowns
- `fixed::Vector{Int}`: global indices of support reactions
- `inactive::Vector{Int}`: excluded slots
- `global_to_free::Vector{Int}`: inverse map — `0` unless free, else the
  position in `free` (i.e. the row/column in the reduced system)
"""
struct DofPartition
    n_global::Int
    free::Vector{Int}
    fixed::Vector{Int}
    inactive::Vector{Int}
    global_to_free::Vector{Int}
end

"first global DOF slot of a node (its Tx); the node owns 6 consecutive slots"
node_dof_start(node::Node) = 6 * (node.index - 1)

"global DOF indices of the element's slots (nodal, then any internal blocks)"
function element_global_dofs(el::AbstractElement)
    g = Vector{Int}(undef, 12)
    s1 = node_dof_start(el.nodeStart)
    s2 = node_dof_start(el.nodeEnd)
    @inbounds for i in 1:6
        g[i] = s1 + i
        g[6+i] = s2 + i
    end
    return g
end

function element_global_dofs(el::VariableElement)
    ni = n_internal_dofs(el)
    @assert el.internal_offset > 0 "VariableElement internal DOFs not yet allocated — call process!"
    g = Vector{Int}(undef, 12 + ni)
    s1 = node_dof_start(el.nodeStart)
    s2 = node_dof_start(el.nodeEnd)
    @inbounds for i in 1:6
        g[i] = s1 + i
        g[6+i] = s2 + i
    end
    @inbounds for i in 1:ni
        g[12+i] = el.internal_offset + i
    end
    return g
end

"""
    partition_dofs(model) -> DofPartition

Classify all global DOFs (see [`DofPartition`](@ref)). Activity is
accumulated from every element's [`dof_signature`](@ref) and every nodal
spring's nonzero stiffness entries; internal element DOFs are always
active. Requires node/element indices to be assigned (done by `process!`).
"""
function partition_dofs(model::Model)
    n_nodal = 6 * length(model.nodes)
    n_internal = sum(n_internal_dofs, model.elements; init=0)
    n_global = n_nodal + n_internal

    active = falses(n_global)
    for el in model.elements
        sig = dof_signature(el)
        g = element_global_dofs(el)
        @inbounds for i in 1:12
            sig[i] && (active[g[i]] = true)
        end
    end
    active[n_nodal+1:end] .= true      # internal DOFs always active

    for sp in model.springs
        s = node_dof_start(sp.node)
        @inbounds for i in 1:6
            sp.stiffness[i] > 0 && (active[s+i] = true)
        end
    end

    free = Int[]
    fixed = Int[]
    inactive = Int[]
    global_to_free = zeros(Int, n_global)
    for (i, node) in enumerate(model.nodes)
        s = 6 * (i - 1)
        for j in 1:6
            g = s + j
            if !node.fixity[j]
                push!(fixed, g)
            elseif active[g]
                push!(free, g)
                global_to_free[g] = length(free)
            else
                push!(inactive, g)
            end
        end
    end
    for g in n_nodal+1:n_global        # internal DOFs are always free
        push!(free, g)
        global_to_free[g] = length(free)
    end

    return DofPartition(n_global, free, fixed, inactive, global_to_free)
end
