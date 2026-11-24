"""
    ModelGeo <: AbstractGeo

Plot-ready arrays extracted from a solved frame `Model`: node positions,
nodal displacements, element connectivity, per-element end-force ranges
(P, Vy, Vz, Tx, My, Mz with their maxima), section properties, and element
vectors. Construct with `Geo(model)` or `ModelGeo(model)`.

Covers the model's frame and truss elements (in model order); surface
elements are skipped — extract them with [`ShellGeo`](@ref).
"""
struct ModelGeo <: AbstractGeo
    nodes::Vector{Vector{Float64}}
    nodes_xy::Vector{Vector{Float64}}
    disp::Vector{Vector{Float64}}
    disp_xy::Vector{Vector{Float64}}
    indices::Vector{Vector{Int64}}
    indices_flat::Vector{Int64}
    P::Vector{Float64}
    max_abs_P::Float64
    Vy::Vector{Float64}
    max_abs_Vy::Float64
    Vz::Vector{Float64}
    max_abs_Vz::Float64
    Tx::Vector{Float64}
    max_abs_Tx::Float64
    My::Vector{Float64}
    max_abs_My::Float64
    Mz::Vector{Float64}
    max_abs_Mz::Float64
    areas::Vector{Float64}
    max_area::Float64
    Ix::Vector{Float64}
    max_Ix::Float64
    Iy::Vector{Float64}
    max_Iy::Float64
    J::Vector{Float64}
    max_J::Float64
    lengths::Vector{Float64}
    element_vectors::Vector{Vector{Float64}}
    element_vectors_xy::Vector{Vector{Float64}}

    function ModelGeo(model::Model)

        results = model.results

        nodes = [Vector(node.position) for node in model.nodes]
        nodes_xy = [node[1:2] for node in nodes]

        disp = [displacement(results, node)[1:3] for node in model.nodes]
        disp_xy = [d[1:2] for d in disp]

        lines = AbstractElement[el for el in model.elements if el isa Union{FrameElement,TrussElement}]
        indices = [[element.nodeStart.index, element.nodeEnd.index] for element in lines]
        indices_flat = isempty(indices) ? Int[] : vcat(indices...)

        forcevectors = [element_forces(results, element) for element in lines]

        P = isempty(forcevectors) ? Float64[] : vcat([forces[[1,7]] .* [-1, 1] for forces in forcevectors]...)
        max_abs_P = maximum(abs, P; init = 0.0)

        Vy = isempty(forcevectors) ? Float64[] : vcat([forces[[2,8]] .* [-1, 1] for forces in forcevectors]...)
        max_abs_Vy = maximum(abs, Vy; init = 0.0)

        Vz = isempty(forcevectors) ? Float64[] : vcat([forces[[3,9]] .* [-1, 1] for forces in forcevectors]...)
        max_abs_Vz = maximum(abs, Vz; init = 0.0)

        Tx = isempty(forcevectors) ? Float64[] : vcat([forces[[4,10]] .* [-1, 1] for forces in forcevectors]...)
        max_abs_Tx = maximum(abs, Tx; init = 0.0)

        My = isempty(forcevectors) ? Float64[] : vcat([forces[[5,11]] .* [-1, 1] for forces in forcevectors]...)
        max_abs_My = maximum(abs, My; init = 0.0)

        Mz = isempty(forcevectors) ? Float64[] : vcat([forces[[6,12]] .* [-1, 1] for forces in forcevectors]...)
        max_abs_Mz = maximum(abs, Mz; init = 0.0)

        sections = [el.section for el in lines]

        areas = Float64[s.A for s in sections]
        max_area = maximum(areas; init = 0.0)

        Ix = Float64[s.Ix for s in sections]
        max_Ix = maximum(Ix; init = 0.0)

        Iy = Float64[s.Iy for s in sections]
        max_Iy = maximum(Iy; init = 0.0)

        J = Float64[s.J for s in sections]
        max_J = maximum(J; init = 0.0)

        element_vectors = [Vector(local_frame(element)[1, :]) for element in lines]
        element_vectors_xy = [evec[1:2] for evec in element_vectors]

        return new(
            nodes,
            nodes_xy,
            disp,
            disp_xy,
            indices,
            indices_flat,
            P,
            max_abs_P,
            Vy,
            max_abs_Vy,
            Vz,
            max_abs_Vz,
            Tx,
            max_abs_Tx,
            My,
            max_abs_My,
            Mz,
            max_abs_Mz,
            areas,
            max_area,
            Ix,
            max_Ix,
            Iy,
            max_Iy,
            J,
            max_J,
            Float64[length(el) for el in lines],
            element_vectors,
            element_vectors_xy
        )
    end

end
