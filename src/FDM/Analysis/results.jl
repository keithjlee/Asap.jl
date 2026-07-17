"""
    NetworkResults{T}

Per-solve results of an FDM network, queried through accessors (the solved
GEOMETRY lives on the nodes — form-finding's deliverable — but forces and
reactions live here):

- `member_force(res, element)` / `res.forces[el.index]` — axial force `q·L`
  [force] at the equilibrium geometry (tension + for positive q)
- `reaction(res, node)` — anchor force on the node's FIXED axes [force]
  (free-axis components are zero: equilibrium holds there)
- `lengths` — member lengths at equilibrium [length]
"""
struct NetworkResults{T}
    lengths::Vector{T}
    forces::Vector{T}
    reactions::Matrix{T}
end

function _network_results(network::Network{T}) where {T}
    cache = network.cache::NetworkCache{T}
    ne = Base.length(network.elements)
    lengths = Vector{T}(undef, ne)
    fs = Vector{T}(undef, ne)
    for (i, el) in enumerate(network.elements)
        lengths[i] = Base.length(el)
        fs[i] = lengths[i] * el.q
    end

    reactions = zeros(T, Base.length(network.nodes), 3)
    rows = rowvals(cache.C)                     # C is ne × nn: column i = node i
    vals = nonzeros(cache.C)
    for (i, node) in enumerate(network.nodes)
        all(node.fixity) && continue
        rxn = zeros(T, 3)
        for p in nzrange(cache.C, i)            # elements at node i
            idx = rows[p]
            el = network.elements[idx]
            eforce = local_x(el) * fs[idx]
            if vals[p] < 0
                rxn .+= eforce
            else
                rxn .-= eforce
            end
        end
        for a in 1:3
            node.fixity[a] && (rxn[a] = zero(T))   # free axis: not a reaction
        end
        reactions[i, :] = rxn
    end

    return NetworkResults{T}(lengths, fs, reactions)
end

"""
    member_force(res::NetworkResults, element::FDMelement) -> T

Axial force in `element` at the solved geometry [force].
"""
member_force(res::NetworkResults, element::FDMelement) = res.forces[element.index]

"""
    reaction(res::NetworkResults, node::FDMnode) -> SVector{3}

Anchor force on `node`'s fixed axes [force]; zero on free axes.
"""
reaction(res::NetworkResults{T}, node::FDMnode) where {T} =
    SVector{3,T}(res.reactions[node.index, :])
