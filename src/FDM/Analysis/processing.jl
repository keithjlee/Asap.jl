"""
    populate_indices!(network::Network)

Populate the node/element global IDs in a network.
"""
function populate_indices!(network::Network)
    for (i, node) in enumerate(network.nodes)
        node.nodeID = i
    end

    for (i, element) in enumerate(network.elements)
        element.elementID = i
        element.iStart = element.pStart.nodeID
        element.iEnd = element.pEnd.nodeID
    end
end


"""
    dofs(points::Vector{FDMnode})

Return the fixed/free DOF information from a vector of nodes.
"""
function dofs(points::Vector{FDMnode})
    getproperty.(points, :dof)
end

"""
Extracts the indices of free (N) and fixed (F) nodes
"""
function NF(d::Union{Vector{Bool}, BitVector})
    return findall(d), findall(.!d)
end

"""
    NF(points::Vector{FDMnode})

Return fixed/free DOF indices from a vector of nodes.
"""
function NF(points::Vector{FDMnode})
    d = dofs(points)
    return NF(d)
end

"""
    NF!(network::Network)

Extract fixed/free DOF indices from a vector of nodes.
"""
function NF!(network::Network)
    network.N, network.F = NF(network.nodes)
end

"""
    branch_matrix(elements::Vector{FDMelement}, points::Vector{FDMnode})

Return the [n_elements × n_nodes] sparse connectivity matrix, `network.C`
"""
function branch_matrix(elements::Vector{FDMelement}, points::Vector{FDMnode})

    #initialize
    c = spzeros(Int64, length(elements), length(points))

    #rows = elements, columns = nodes
    for (i, element) in enumerate(elements)
        c[i, element.iStart] = -1
        c[i, element.iEnd] = 1
    end

    return c
end

"""
    branch_matrix!(network::Network)

Populate the [n_elements × n_nodes] sparse connectivity matrix, `network.C`
"""
function branch_matrix!(network::Network)
    network.C = branch_matrix(network.elements, network.nodes)
end

"""
    deconstruct_nodes(points::Vector{FDMnode})

return the [n_nodes × 3] matrix of nodal positions, `network.xyz`
"""
function deconstruct_nodes(points::Vector{FDMnode})
    positions = getproperty.(points, :position)
    return Float64.([getindex.(positions, 1) getindex.(positions, 2) getindex.(positions, 3)])
end

"""
    deconstruct_nodes!(network::Network)

Populate the [n_nodes × 3] matrix of nodal positions, `network.xyz`
"""
function deconstruct_nodes!(network::Network)
    network.xyz = deconstruct_nodes(network.nodes)
end

"""
    force_densities(elements::Vector{FDMelement})

Extract the force densities of a vector of elements.
"""
function force_densities(elements::Vector{FDMelement})
    getproperty.(elements, :q)
end

"""
    force_densities!(network::Network)

Extract the force densities of elements in a network and assemble the global Q matrix, `network.Q`
"""
function force_densities!(network::Network)
    network.q = force_densities(network.elements)
    network.Q = spdiagm(network.q)
end

"""
    load_matrix!(network::Network)

Generate the [n_nodes × 3] matrix of nodal forces, P
"""
function load_matrix!(network::Network)
    network.P = zeros(length(network.nodes), 3)

    for load in network.loads
        network.P[load.point.nodeID, :] = load.force
    end
    
end


"""
    process!(network::Network)

Preprocess a network. Assign indices and references, and generate global collectors of position, force density, etc.
"""
function process!(network::Network)
    #populate indices
    populate_indices!(network)

    #fixed-free indices N, F
    NF!(network)

    #branch node matrix C
    branch_matrix!(network)
    network.Cn = network.C[:, network.N]
    network.Cf = network.C[:, network.F]

    #nodal positions [x y z]
    deconstruct_nodes!(network)

    #force density vector q
    force_densities!(network)

    #load matrix P
    load_matrix!(network)
    network.Pn = network.P[network.N, :]

    network.processed = true
end