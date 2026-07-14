"""
    TrussGeo

Utility structure for visualizing truss models
"""
struct TrussGeo <: AbstractGeo
    nodes::Vector{Vector{Float64}}
    nodes_xy::Vector{Vector{Float64}}
    disp::Vector{Vector{Float64}}
    disp_xy::Vector{Vector{Float64}}
    indices::Vector{Vector{Int64}}
    indices_flat::Vector{Int64}
    forces::Vector{Float64}
    max_abs_force::Float64
    areas::Vector{Float64}
    max_area::Float64
    lengths::Vector{Float64}
    element_vectors::Vector{Vector{Float64}}
    element_vectors_xy::Vector{Vector{Float64}}
    load_positions::Vector{Vector{Float64}}
    load_positions_xy::Vector{Vector{Float64}}
    load_vectors::Vector{Vector{Float64}}
    load_vectors_xy::Vector{Vector{Float64}}

    function TrussGeo(model::Model)

        results = model.results

        nodes = [Vector(node.position) for node in model.nodes]
        nodes_xy = [node[1:2] for node in nodes]

        disp = [displacement(results, node)[1:3] for node in model.nodes]
        disp_xy = [d[1:2] for d in disp]

        indices = [[element.nodeStart.index, element.nodeEnd.index] for element in model.elements]
        indices_flat = vcat(indices...)

        forces = [axial_force(results, element) for element in model.elements]
        max_abs_force = maximum(abs.(forces))

        areas = getproperty.(getproperty.(model.elements, :section), :A)
        max_area = maximum(areas)

        load_positions = [Vector(load.node.position) for load in model.loads]
        load_positions_xy = [load[1:2] for load in load_positions]

        load_vectors = [Vector(load.value) for load in model.loads]
        load_vectors_xy = [load[1:2] for load in load_vectors]

        element_vectors = [Vector(local_frame(element)[1, :]) for element in model.elements]
        element_vectors_xy = [evec[1:2] for evec in element_vectors]


        return new(
            nodes,
            nodes_xy,
            disp,
            disp_xy,
            indices,
            indices_flat,
            forces,
            max_abs_force,
            areas,
            max_area,
            length.(model.elements),
            element_vectors,
            element_vectors_xy,
            load_positions,
            load_positions_xy,
            load_vectors,
            load_vectors_xy
        )

    end

end
