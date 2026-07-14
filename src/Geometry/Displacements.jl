#=
Displaced-shape sampling for plotting.

v1.0: implemented as a thin layer over Asap's exact displacement recovery —
`internal_forces(model, element)` builds a closed-form element state and
`local_displacements(state, t)` evaluates the exact local (u, v, w) at any
fraction t. This replaces the legacy shape-function + per-load catalog
implementation.
=#

"""
    ElementDisplacements

The displaced shape of one element, sampled at evenly spaced stations using
the exact elastic curve from Asap's displacement recovery.

# Fields
- `element::AbstractElement` the sampled element
- `resolution::Integer` number of stations
- `x::Vector` station positions along the member [length]
- `ulocal::Matrix` [3 × n] displacements in the element LCS (axial, local y, local z) [length]
- `uglobal::Matrix` [3 × n] displacements in the GCS [length]
- `basepositions::Matrix` [3 × n] undisplaced sampling positions [length]
"""
struct ElementDisplacements
    element::AbstractElement
    resolution::Integer
    x::Vector{Float64}
    ulocal::Matrix{Float64}
    uglobal::Matrix{Float64}
    basepositions::Matrix{Float64}
end

"""
    ElementDisplacements(element::FrameElement, model::Model; resolution = 20)

Get the local/global displacements of an element.

- `ulocal` [3 × resolution] displacements in the element LCS (axial, local y, local z)
- `uglobal` [3 × resolution] displacements in the GCS
- `basepositions` [3 × resolution] undisplaced sampling positions along the element
"""
function ElementDisplacements(element::FrameElement, model::Model; resolution = 20)

    L = length(element)
    tinc = collect(range(0, 1, resolution))
    xinc = tinc .* L

    Λ = local_frame(element)
    state = internal_forces(model, element)

    ulocal = zeros(3, resolution)
    for (i, t) in enumerate(tinc)
        ulocal[:, i] = local_displacements(state, t)
    end

    # rows of Λ are the local axes in GCS: global = Λ' × local
    uglobal = Matrix(Λ' * ulocal)

    basepositions = Vector(element.nodeStart.position) .+ Vector(Λ[1, :]) * xinc'

    return ElementDisplacements(element, resolution, xinc, ulocal, uglobal, basepositions)
end

"""
    ElementDisplacements(element::TrussElement, model::Model; resolution = 20)

Truss elements displace linearly between their end nodes.
"""
function ElementDisplacements(element::TrussElement, model::Model; resolution = 20)

    L = length(element)
    tinc = collect(range(0, 1, resolution))
    xinc = tinc .* L

    Λ = local_frame(element)

    u1 = displacement(model.results, element.nodeStart)[1:3]
    u2 = displacement(model.results, element.nodeEnd)[1:3]

    uglobal = hcat([(1 - t) .* u1 .+ t .* u2 for t in tinc]...)
    ulocal = Matrix(Λ * uglobal)

    basepositions = Vector(element.nodeStart.position) .+ Vector(Λ[1, :]) * xinc'

    return ElementDisplacements(element, resolution, xinc, ulocal, uglobal, basepositions)
end

"""
    ElementDisplacements(elements::Vector{<:AbstractElement}, model::Model; resolution = 20)

Get the displacements of an ordered chain of elements that form a single physical member.
"""
function ElementDisplacements(elements::Vector{<:AbstractElement}, model::Model; resolution = 20)

    xstore = Vector{Float64}()
    ulocalstore = Vector{Matrix{Float64}}()
    uglobalstore = Vector{Matrix{Float64}}()
    basepointstore = Vector{Matrix{Float64}}()

    resolution = max(Int(round(resolution / length(elements))), 2)

    for element in elements

        ed = ElementDisplacements(element, model; resolution = resolution)

        if isempty(xstore)
            xstore = [xstore; ed.x]
        else
            xstore = [xstore; xstore[end] .+ ed.x]
        end

        push!(ulocalstore, ed.ulocal)
        push!(uglobalstore, ed.uglobal)
        push!(basepointstore, ed.basepositions)
    end

    return ElementDisplacements(elements[1], resolution, xstore, hcat(ulocalstore...), hcat(uglobalstore...), hcat(basepointstore...))
end

"""
    displacements(model::Model, increment::Real)

Get the displacements of all elements in a model
"""
function displacements(model::Model, increment::Real)
    results = Vector{ElementDisplacements}()

    for element in model.elements
        L = length(element)
        n = max(Int(round(L / increment)), 2)

        push!(results, ElementDisplacements(element, model; resolution = n))
    end

    return results
end
