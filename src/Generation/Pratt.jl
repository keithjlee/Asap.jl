"""
    Pratt2D <: AbstractGenerator

A 2D Pratt truss in the XY plane, generated and solved on construction.

# Fields
- `model::Model{Float64}` solved structural model
- `L::Real` total span [length]
- `nbays::Integer` number of bays
- `dx::Real` bay spacing = L / nbays [length]
- `dy::Real` truss depth [length]
"""
struct Pratt2D <: AbstractGenerator
    model::Model{Float64}
    L::Real
    nbays::Integer
    dx::Real
    dy::Real
end

function Pratt2D(
    L::Real,
    nbays::Integer,
    dy::Real,
    section::AbstractSection;
    load = [0., -1., 0.],
    up = true)

    #make sure nbays is even
    @assert nbays % 2 == 0 "nbays must be even for a pratt truss"

    #x positions of nodes
    xrange = range(0, L, nbays+1)

    #make main nodes
    main_chord_nodes = [Node([x, 0., 0.], :free, :main) for x in xrange]
    n_main = nbays + 1
    
    #make supports
    fixnode!(first(main_chord_nodes), :pinned)
    first(main_chord_nodes).id = :pin

    fixnode!(last(main_chord_nodes), :xfree)
    last(main_chord_nodes).id = :roller

    #make offset nodes
    offset = up ? dy : -dy
    offset_chord_nodes = [Node([x, offset, 0.], :free, :offset) for x in xrange[2:end-1]]
    n_offset = length(offset_chord_nodes)

    #make main chord elements
    main_chord_elements = [TrussElement(main_chord_nodes[i], main_chord_nodes[i+1], section, :main) for i = 1:nbays]

    #main offset chord elements
    offset_chord_elements = [TrussElement(offset_chord_nodes[i], offset_chord_nodes[i+1], section, :offset) for i = 1:n_offset-1]

    #make vertical strut elements
    strut_elements = [TrussElement(n1, n2, section, :strut) for (n1, n2) in zip(main_chord_nodes[2:end-1], offset_chord_nodes)]

    #make web elements
    i_mid_main = div(n_main, 2, RoundUp)
    i_mid_offset = div(n_offset, 2, RoundUp)

    #left side
    i_left_main = 1:i_mid_main-1
    i_left_offset = 1:i_mid_offset

    web_left = [TrussElement(main_chord_nodes[i1], offset_chord_nodes[i2], section, :web) for (i1, i2) in zip(i_left_main, i_left_offset)]

    #right side
    i_right_main = n_main:-1:i_mid_main+1
    i_right_offset = n_offset:-1:i_mid_offset

    web_right = [TrussElement(main_chord_nodes[i1], offset_chord_nodes[i2], section, :web) for (i1, i2) in zip(i_right_main, i_right_offset)]

    #collect
    nodes = [main_chord_nodes; offset_chord_nodes]
    elements = [main_chord_elements; offset_chord_elements; strut_elements; web_left; web_right]

    #loads
    loads = [NodeForce(node, load) for node in nodes[:main]]

    #assemble
    model = Model(nodes, elements, loads)
    planarize!(model)
    solve!(model)

    return Pratt2D(model, L, nbays, xrange[2], dy)

end