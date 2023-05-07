# set_theme!(kjl_light)

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
Creates a 1x3 vector of nodal position
"""
function deconstructnode(point::FDMnode; horizontal = true)
    if horizontal
        return Float64.([point.x point.y point.z])
    else
        return Float64.([point.x, point.y, point.z])
    end
end

"""
Creates an nx3 matrix of nodal positions
"""
function deconstructnode(points::Vector{FDMnode})
    return Float64.([getproperty.(points, :x) getproperty.(points, :y) getproperty.(points, :z)])
end

"""
Directly on a network
"""
function deconstructnodes!(network::Network)
    network.xyz = deconstructnode(network.nodes)
end


"""
Creates load matrix P (nx3) with loads indexed vertically w/r/t node order
"""
function loadMatrix(loads::Vector{FDMload}, n::Int64)
    p = zeros(n, 3)
    for load in loads
        p[load.index, :] = load.force
    end
    return p
end

"""
directly on a network
"""
function loadMatrix!(network::Network)
    network.P = loadMatrix(network.loads, length(network.nodes))
end



"""
Preprocessing of data for a new FDM network
"""
function process!(network::Network)
    #fixed-free indices N, F
    NF!(network)

    #branch node matrix C
    branchmatrix!(network)
    network.Cn = network.C[:, network.N]
    network.Cf = network.C[:, network.F]

    #nodal positions [x y z]
    deconstructnodes!(network)

    #force density vector q
    forcedensities!(network)

    #load matrix P
    loadMatrix!(network)
    network.Pn = network.P[network.N, :]

    network.processed = true
end

"""
performs the fdm analysis testing
"""
function solve!(network::Network; reprocess = false)
    if !network.processed || reprocess
        process!(network)
    end

    # solve for final free node positions
    network.xyz[network.N, :] = (network.Cn' * network.Q * network.Cn) \ (network.Pn - network.Cn' * network.Q * network.Cf * network.xyz[network.F, :])

    # update nodal positions
    xyzupdate!(network)
end

"""
Solve w/r/t q ...
"""
function solve(network::Network, q::Union{Vector{Int64}, Vector{Float64}})

    @assert length(q) == length(network.elements) "q and elements must have equal length"

    Q = spdiagm(q)

    xyzout = copy(network.xyz)
    xyzout[network.N, :] .= (network.Cn' * Q * network.Cn) \ (network.Pn - network.Cn' * Q * network.Cf * network.xyz[network.F, :])

    return xyzout
end

"""
Extracts force density vector (q) from elements
"""
function forcedensities(elements::Vector{FDMelement})
    getproperty.(elements, :q)
end

"""
Extracts q from elements and populates network fields
"""
function forcedensities!(network::Network)
    network.q = forcedensities(network.elements)
    network.Q = spdiagm(network.q)
end

"""
Creates the branch-node connectivity matrix C
"""
function branchmatrix(elements::Vector{FDMelement}, points::Vector{FDMnode})

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
function branchmatrix!(network::Network)
    network.C = branchmatrix(network.elements, network.nodes)
end

"""
∑|F|L
"""
function FL(network::Network)
    return sum(norm.(eachrow(network.C * network.xyz)).^2 .* network.q)
end

"""
update nodal positions
"""
function vec2node!(vec::Vector{<:Real}, node::FDMnode)
    node.x, node.y, node.z = vec
end

"""
applied to network
"""
function xyzupdate!(network::Network)
    for index in network.N
        vec2node!(network.xyz[index, :], network.nodes[index])
    end
end