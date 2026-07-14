"""
    Warren2D <: AbstractGenerator

A 2D Warren truss (arch or catenary form) in the XY plane, generated and
solved on construction.

# Fields
- `model::Model{Float64}` solved structural model
- `n::Integer` number of bays (odd)
- `dx::Real` bay span [length]
- `dy::Real` truss depth [length]
- `section::AbstractSection` cross section applied to all members
- `type::Symbol` `:arch` (long chord at bottom) or `:catenary` (long chord at top)
"""
struct Warren2D <: AbstractGenerator
    model::Model{Float64}
    n::Integer
    dx::Real
    dy::Real
    section::AbstractSection
    type::Symbol
end

"""
    generatewarren2d(n::Integer,...)

Generate a 2D warren truss in the XY plane.

Required inputs:
- `n::Integer` Number of bays
- `dx::Real` Bay span
- `dy::Real` truss depth
- `section::AbstractSection` cross section of elements

Default inputs:
- `type::Symbol = :arch` :arch = long chord at bottom; :catenary = long chord at top
"""
function Warren2D(n::Integer,
        dx::Real,
        dy::Real,
        section::AbstractSection;
        load = [0., -1., 0.],
        type = :arch)

    @assert n % 2 != 0 "n must be odd"
    @assert type == :arch || type == :catenary "type must be :arch or :catenary"

    #counters
    count = 1
    longids = Vector{Int64}()

    #node collector
    nodes = Vector{Node{Float64}}()

    #generate longer chord first
    if type == :arch
        longid = :bottomchord
        shortid = :topchord
        y = dy
    else
        longid = :topchord
        shortid = :bottomchord
        y = -dy
    end

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
            node.id = longid
        end

        push!(nodes, node)
        push!(longids, count)

        count += 1
    end

    #generate shorter chord
    shortids = Vector{Int64}()
    x0 = dx / 2

    for i = 1:n-1
        xposition = x0 + dx * (i - 1)

        node = Node([xposition, y, 0.], :free)
        node.id = shortid

        push!(nodes, node)
        push!(shortids, count)
        count += 1
    end

    #elements
    elements = Vector{AbstractElement{Float64}}()

    #long chords
    for i = 1:n-1
        element = TrussElement(nodes[longids[i:i+1]]..., section)
        element.id = longid

        push!(elements, element)
    end

    #short chords
    for i = 1:n-2
        element = TrussElement(nodes[shortids[i:i+1]]..., section)
        element.id = shortid

        push!(elements, element)
    end

    #webs
    for i = 1:n-1
        element = TrussElement(nodes[[longids[i], shortids[i]]]..., section)
        element.id = :web
        push!(elements, element)

        element = TrussElement(nodes[[shortids[i], longids[i+1]]]..., section)
        element.id = :web
        push!(elements, element)
    end

    #dummy load
    loads = [NodeForce(n, load) for n in nodes[longid]]

    #assemble and solve
    model = Model(nodes, elements, loads)
    planarize!(model)
    solve!(model)

    #collect data
    truss = Warren2D(model, n, dx, dy, section, type)

    #output
    return truss
end

function Warren2D(xpositions::Vector{<:Real},
    ypositions::Vector{<:Real},
    ypositions2::Vector{<:Real},
    section::AbstractSection;
    type = :arch)

    @assert length(xpositions) == length(ypositions) == length(ypositions2) + 1
    @assert type == :arch || type == :catenary "type must be :arch or :catenary"

    #counters
    count = 1
    longids = Vector{Int64}()

    #node collector
    nodes = Vector{Node{Float64}}()

    #generate longer chord first
    if type == :arch
        longid = :bottomchord
        shortid = :topchord
    else
        longid = :topchord
        shortid = :bottomchord
    end

    n = length(xpositions)
    i = 1

    ## generate long chord up to symmetry
    for (x, y) in zip(xpositions, ypositions)

        node = Node([x, y, 0.], :free)
        if i == 1
            fixnode!(node, :pinned)
            node.id = :pin
            i += 1
        else
            node.id = longid
        end

        push!(nodes, node)
        push!(longids, count)
        count += 1
    end

    ## generate other side of symmetry
    Lhalf = last(xpositions)
    incs = Lhalf .- reverse(xpositions[1:end-1])
    for (inc, y) in zip(incs, reverse(ypositions[1:end-1]))
        node = Node([inc, y, 0.], :free)
        node.id = longid
        push!(nodes, node)
        push!(longids, count)
        count += 1
    end
    fixnode!(last(nodes), :yfixed)
    last(nodes).id = :roller


    #generate shorter chord
    shortids = Vector{Int64}()
    for i = 1:length(ypositions2)

        xposition = (xpositions[i] + xpositions[i+1]) / 2

        node = Node([xposition, ypositions2[i], 0.], :free)
        node.id = shortid

        push!(nodes, node)
        push!(shortids, count)
        count += 1
    end

    #other side
    for (i, y) in zip(length(xpositions):length(longids), reverse(ypositions2))
        x = (nodes[i].position[1] + nodes[i+1].position[1]) / 2

        node = Node([x, y, 0.], :free)
        node.id = shortid

        push!(nodes, node)
        push!(shortids, count)

        count += 1
    end


    #elements
    elements = Vector{AbstractElement{Float64}}()

    #long chords
    for i = 1:length(longids) - 1
        element = TrussElement(nodes[longids[i:i+1]]..., section)
        element.id = longid

        push!(elements, element)
    end

    #short chords
    for i = 1:length(shortids) - 1
        element = TrussElement(nodes[shortids[i:i+1]]..., section)
        element.id = shortid

        push!(elements, element)
    end

    #webs
    for i = 1:length(longids) - 1
        element = TrussElement(nodes[[longids[i], shortids[i]]]..., section)
        element.id = :web
        push!(elements, element)

        element = TrussElement(nodes[[shortids[i], longids[i+1]]]..., section)
        element.id = :web
        push!(elements, element)
    end

    #dummy load
    loads = [NodeForce(n, [0., -1., 0.],) for n in nodes[longid]]

    #assemble and solve
    model = Model(nodes, elements, loads)
    planarize!(model)
    solve!(model)

    #collect data
    #NOTE: pre-existing breakage (unported): `dx` and `dy` were never defined in this
    #method prior to the v1.0 port; this constructor has always thrown at this line.
    truss = Warren2D(model, n, dx, dy, section, type)

end
