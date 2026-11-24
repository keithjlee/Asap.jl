"""
    clear_supports!(model::Model)

Remove all fixed DOFs from a model and reprocess. Will require new definitions of fixed supports before running `solve!()`
"""
function clear_supports!(model::Model)

    for node in model.nodes
        fixnode!(node, :free)
        node.id = :free
    end

    process!(model)

end

"""
    element_connectivity(model::Model)

Get the Cel = [n_element × n_element] adjacency matrix for a given model, where Cel[i, j] = 1 if element i shares a node with element j, and is 0 otherwise.
"""
function element_connectivity(model::Model)
    n_elements = length(model.elements)
    # node → elements touching it (any element type: nodes(el))
    touching = Dict{Int,Vector{Int}}()
    for (e, el) in enumerate(model.elements), n in nodes(el)
        push!(get!(touching, n.index, Int[]), e)
    end
    I_ = Int[]; J_ = Int[]
    for (_, els) in touching, a in els, b in els
        a == b && continue
        push!(I_, a); push!(J_, b)
    end
    Cel = sparse(I_, J_, ones(Int64, length(I_)), n_elements, n_elements)
    return min.(Cel, 1)   # an element pair sharing two nodes still counts once
end
