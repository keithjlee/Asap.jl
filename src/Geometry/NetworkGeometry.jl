"""
    NetworkGeo <: AbstractGeo

Plot-ready arrays extracted from a solved FDM `Network`: node positions,
element connectivity, force densities and member forces/lengths with their
maxima. Construct with `Geo(network)` or `NetworkGeo(network)`.
"""
struct NetworkGeo <: AbstractGeo
    nodes::Vector{Vector{Float64}}
    indices::Vector{Vector{Int64}}
    indices_flat::Vector{Int64}
    q::Vector{Float64}
    max_abs_q::Float64
    forces::Vector{Float64}
    max_abs_force::Float64
    lengths::Vector{Float64}
    element_vectors::Vector{Vector{Float64}}
    load_positions::Vector{Vector{Float64}}
    load_vectors::Vector{Vector{Float64}}

    function NetworkGeo(network::Network)
        nodes = [Vector(row) for row in eachrow(network.xyz)]
        indices = [[element.iStart, element.iEnd] for element in network.elements]
        indices_flat = vcat(indices...)

        q = network.q
        max_abs_q = maximum(abs.(q))

        forces = force.(network.elements)
        max_abs_force = maximum(abs.(forces))

        lengths = length.(network.elements)

        element_vectors = local_x.(network.elements; unit = false)

        load_positions = [load.point.position for load in network.loads]
        load_vectors = getproperty.(network.loads, :force)

        new(
            nodes,
            indices,
            indices_flat,
            q,
            max_abs_q,
            forces,
            max_abs_force,
            lengths,
            element_vectors,
            load_positions,
            load_vectors
        )

    end

end