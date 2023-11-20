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
Extracts the vector of Boolean DOFs in node sequence.
|dofs| = |nodes|
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
Extracts the fixed/free indices of a vector of nodes
"""
function NF(points::Vector{FDMnode})
    d = dofs(points)
    return NF(d)
end

"""
Directly on network
"""
function NF!(network::Network)
    network.N, network.F = NF(network.nodes)
end

"""
Creates the branch-node connectivity matrix C
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
Directly on a network
"""
function branch_matrix!(network::Network)
    network.C = branch_matrix(network.elements, network.nodes)
end

"""
Creates an nx3 matrix of nodal positions
"""
function deconstruct_nodes(points::Vector{FDMnode})
    positions = getproperty.(points, :position)
    return Float64.([getindex.(positions, 1) getindex.(positions, 2) getindex.(positions, 3)])
end

"""
Directly on a network
"""
function deconstruct_nodes!(network::Network)
    network.xyz = deconstruct_nodes(network.nodes)
end

"""
Extracts force density vector (q) from elements
"""
function force_densities(elements::Vector{FDMelement})
    getproperty.(elements, :q)
end

"""
Extracts q from elements and populates network fields
"""
function force_densities!(network::Network)
    network.q = force_densities(network.elements)
    network.Q = spdiagm(network.q)
end

"""
directly on a network
"""
function load_matrix!(network::Network)
    network.P = zeros(length(network.nodes), 3)

    for load in network.loads
        network.P[load.point.nodeID, :] = load.force
    end
    
end


"""
Preprocessing of data for a new FDM network
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