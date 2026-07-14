"""
    BakerTruss <: AbstractGenerator

A 2D modified Baker truss with X-bracing in each bay offset to quarter points.

# Fields
- `model::Model{Float64}` structural model
- `nbays::Integer` number of bays across span
- `dx::Real` x-spacing of bays
- `dy::Real` y-spacing of bays
- `section::AbstractSection` section applied to elements
"""
struct BakerTruss <: AbstractGenerator
    model::Model{Float64}
    nbays::Integer
    dx::Real
    dy::Real
    section::AbstractSection
end

"""
    BakerTruss(nbays::Integer, dx::Real, dy::Real, section::AbstractSection, load = [0., -10., 0.])

Generate a 2D modified Baker truss.

# Arguments
- `nbays::Integer` number of total bays (must be even)
- `dx::Real` x-spacing of bays
- `dy::Real` y-spacing of bays
- `section::AbstractSection` cross section of bars
- `load::Vector{Float64} = [0., -10., 0.]` Load applied to free bottom nodes
"""
function BakerTruss(nbays::Integer, dx::Real, dy::Real, section::AbstractSection, load = [0., -10., 0.])

    @assert nbays % 2 == 0 "n must be even"
    n = nbays + 1

    #generate bottom nodes
    bottomnodes = Vector{Node{Float64}}()
    for i = 1:n
        xposition = dx * (i - 1)

        node = Node([xposition, 0., 0.], :free)
        if i == 1
            fixnode!(node, :pinned)
            node.id = :pin
        elseif i == n
            node.fixity = vcat([true, false, false], trues(3))
            node.id = :roller
        else
            node.id = :bottom
        end

        push!(bottomnodes, node)
    end

    #generate top nodes
    topnodes = Vector{Node{Float64}}()
    for i = 1:n
        xposition = dx * (i - 1)

        node = Node([xposition, dy, 0.], :free, :top)

        push!(topnodes, node)
    end

    #generate web nodes

    #left side
    x0 = .25dx #initial offset

    webleft = Vector{Node{Float64}}()
    for i = 1:nbays/2
        xposition = x0 + (i - 1) * dx

        node = Node([xposition, dy/2, 0.], :free, :web)

        push!(webleft, node)
    end

    #right side
    x0 = dx * (nbays/2 + 0.75)
    webright = Vector{Node{Float64}}()
    for i = 1:nbays/2
        xposition = x0 + (i - 1) * dx

        node = Node([xposition, dy/2, 0.], :free, :web)

        push!(webright, node)
    end

    webnodes = [webleft; webright]

    #elements

    #bottom chord
    bottomchord = [TrussElement(bottomnodes[i], bottomnodes[i+1], section, :bottom) for i = 1:nbays]

    #topchord
    topchord = [TrussElement(topnodes[i], topnodes[i+1], section, :top) for i = 1:nbays]

    #verticals
    verts = [TrussElement(n1, n2, section, :vertical) for (n1, n2) in zip(bottomnodes, topnodes)]

    #webs
    webs = Vector{AbstractElement{Float64}}()
    for i = 1:nbays
        #topleft
        push!(webs, TrussElement(webnodes[i], topnodes[i], section, :web))
        #bottomleft
        push!(webs, TrussElement(webnodes[i], bottomnodes[i], section, :web))

        #topright
        push!(webs, TrussElement(webnodes[i], topnodes[i+1], section, :web))
        #bottomright
        push!(webs, TrussElement(webnodes[i], bottomnodes[i+1], section, :web))
    end

    nodes = [bottomnodes; topnodes; webnodes]
    elements = [bottomchord; topchord; verts; webs]
    loads = [NodeForce(node, load) for node in nodes if node.id == :bottom]

    model = Model(nodes, elements, loads)
    planarize!(model)
    solve!(model)

    BakerTruss(
        model,
        nbays,
        dx,
        dy,
        section
    )
end