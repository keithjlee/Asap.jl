"""
    SpaceFrame <: AbstractGenerator

A double-layer spaceframe grid (offset square-on-square), generated and
solved on construction.

# Fields
- `model::Model{Float64}` solved structural model
- `nx::Integer`, `dx::Real` bays and bay span in x [length]
- `ny::Integer`, `dy::Real` bays and bay span in y [length]
- `dz::Real` layer offset (roof thickness) [length]
- `section::AbstractSection` cross section applied to all members
- `support::Symbol` support layout: `:corner`, `:x`, `:y`, or `:xy`
- `load::Vector` nodal load applied to top-layer nodes [force]
- `base::Vector` origin offset of the generated grid [length]
- `ibottom`, `itop::Matrix{Int64}` node indices of the bottom/top layers
- `isquares::Matrix{Vector{Int64}}` node indices of each top-layer module
- `isupport::Vector{Int64}` supported node indices
- `iX1`, `iX2`, `iY1`, `iY2::Vector{Int64}` node indices of the four exterior edges
"""
struct SpaceFrame <: AbstractGenerator
    model::Model{Float64}
    nx::Integer
    dx::Real
    ny::Integer
    dy::Real
    dz::Real
    section::AbstractSection
    support::Symbol
    load::Vector{<:Real}
    base::Vector{<:Real}
    ibottom::Matrix{Int64}
    itop::Matrix{Int64}
    isquares::Matrix{Vector{Int64}}
    isupport::Vector{Int64}
    iX1::Vector{Int64}
    iX2::Vector{Int64}
    iY1::Vector{Int64}
    iY2::Vector{Int64}
end

"""
    SpaceFrame(nx::Integer,
    dx::Real,
    ny::Integer,
    dy::Real,
    dz::Real,
    section::AbstractSection;
    support = :corner,
    load = [0., 0., -10.],
    base = [0., 0., 0.])

Generate a spaceframe truss model of planar bay spacing dx, dy, with nx, ny nodes in each direction, with roof thickness dz. Options for support include:
- `:corner` pinned supports at the 4 corner bays (total of 16 pinned supports)
- `:x` pinned supports along exterior nodes parallel to the global X direction
- `:y` pinned supports along exterior nodes parallel to the global Y direction
- `:xy` pinned supports along all exterior nodes
"""
function SpaceFrame(nx::Integer,
        dx::Real,
        ny::Integer,
        dy::Real,
        dz::Real,
        section::AbstractSection;
        support = :corner,
        load = [0., 0., -10.],
        base = [0., 0., 0.])

    #generate nodes for bottom plane
    bottomnodes = [Node([dx * (i-1), dy * (j-1), 0.] .+ base, :free) for i in 1:nx+1, j in 1:ny+1]
    for node in bottomnodes
        node.id = :bottom
    end

    #generate top nodes
    xinit = dx / 2
    yinit = dy / 2

    topnodes = [Node([dx * (i-1) + xinit, dy * (j-1) + yinit, dz], :free) for i in 1:nx, j in 1:ny]
    for node in topnodes
        node.id = :top
    end

    #elements
    elements = Vector{AbstractElement{Float64}}()

    #generate bottom horizontal elements
    #parallel to x
    for j = 1:ny+1
        for i = 1:nx
            element = TrussElement(bottomnodes[i,j], bottomnodes[i+1,j], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #parallel to y
    for i = 1:nx+1
        for j = 1:ny
            element = TrussElement(bottomnodes[i,j], bottomnodes[i,j+1], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #generate top horizontal elements
    #parallel to x
    for j = 1:ny
        for i = 1:nx-1
            element = TrussElement(topnodes[i,j], topnodes[i+1,j], section)
            element.id = :top

            push!(elements, element)
        end
    end

    #parallel to y
    for i = 1:nx
        for j = 1:ny-1
            element = TrussElement(topnodes[i,j], topnodes[i,j+1], section)
            element.id = :top

            push!(elements, element)
        end
    end

    #generate web elements
    for i = 1:nx
        for j = 1:ny
            e1 = TrussElement(topnodes[i,j], bottomnodes[i,j], section)
            e2 = TrussElement(topnodes[i,j], bottomnodes[i+1,j], section)
            e3 = TrussElement(topnodes[i,j], bottomnodes[i,j+1], section)
            e4 = TrussElement(topnodes[i,j], bottomnodes[i+1,j+1], section)

            e1.id = e2.id = e3.id = e4.id = :web

            push!(elements, e1, e2, e3, e4)
        end
    end


    #generate node index matrices
    count = 1

    ibottomnodes = zeros(Int64, nx+1, ny+1)
    for j in 1:ny+1
        for i in 1:nx+1
            ibottomnodes[i,j] = count
            count += 1
        end
    end

    itopnodes = zeros(Int64, nx, ny)
    for j in 1:ny
        for i in 1:nx
            itopnodes[i,j] = count
            count += 1
        end
    end

    isquares = [[ibottomnodes[i,j],
        ibottomnodes[i+1,j],
        ibottomnodes[i,j+1],
        ibottomnodes[i+1,j+1]] for i in 1:nx, j in 1:ny]


    #generate supports
    ix1 = ibottomnodes[1,:]
    ix2 = ibottomnodes[end,:]
    iy1 = ibottomnodes[:,1]
    iy2 = ibottomnodes[:,end]

    #generate supports
    if support == :corner
        for i in [1, nx]
            for j in [1, ny]
                idset = isquares[i,j]

                for k in idset
                    fixnode!(bottomnodes[k], :pinned)
                    bottomnodes[k].id = :support
                end
            end
        end
    elseif support == :center
        isupport = nx % 2 == 1 ? [Int(ceil(nx / 2))] : Int.([ceil(nx / 2), floor(nx /2)])
        jsupport = ny % 2 == 1 ? [Int(ceil(ny / 2))] : Int.([ceil(ny / 2), floor(ny /2)])

        iset = vcat([vec(isquares[i,j]) for i in isupport, j in jsupport]...)

        for i in iset
            bottomnodes[i].id = :support
            fixnode!(bottomnodes[i], :pinned)
        end
    elseif support == :x
        for i in [ix1; ix2]
            bottomnodes[i].id = :support
            fixnode!(bottomnodes[i], :pinned)
        end

    elseif support == :y
        for i in [iy1; iy2]
            bottomnodes[i].id = :support
            fixnode!(bottomnodes[i], :pinned)
        end

    elseif support == :xy
        for i in [ix1; ix2; iy1; iy2]
            bottomnodes[i].id = :support
            fixnode!(bottomnodes[i], :pinned)
        end
    else
        error("support must be: :corner, :center, :x, :y, :xy")
    end

    flatnodes = [vec(bottomnodes); vec(topnodes)]

    #generate loads
    loads  = [NodeForce(node, load) for node in flatnodes if node.id != :support]

    #assemble
    model = Model(flatnodes, elements, loads)
    solve!(model)

    isupport = findall(model.nodes, :support)


    spaceframe = SpaceFrame(model,
        nx,
        dx,
        ny,
        dy,
        dz,
        section,
        support,
        load,
        base,
        ibottomnodes,
        itopnodes,
        isquares,
        isupport,
        ix1,
        ix2,
        iy1,
        iy2)

    return spaceframe
end

"""
    SpaceFrame(nx, dx, ny, dy, z0, interpolator, section, offset = false; ...)

Variable-depth spaceframe: the top-layer height is sampled from
`interpolator(u, v)`, any callable mapping normalized plan coordinates
`u, v ∈ [0, 1]` to a height offset added to `z0`. Pass e.g. a closure over
an analytic surface, or an `Interpolations.jl` extrapolation built on a
`(0,1) × (0,1)` grid — the generator only ever calls it as a function.
"""
function SpaceFrame(nx::Integer,
    dx::Real,
    ny::Integer,
    dy::Real,
    z0::Real,
    interpolator,
    section::AbstractSection,
    offset = false;
    support = :corner,
    load = [0., 0., -10.],
    base = [0., 0., 0.])

    xmax = dx * nx
    ymax = dy * ny

    if offset
        #generate nodes for bottom plane
        bottomnodes = [Node([dx * (i-1), 
            dy * (j-1), 
            interpolator(dx * (i-1) / xmax, dy * (j-1) /ymax)] .+ base, :free) for i in 1:nx+1, j in 1:ny+1]
        for node in bottomnodes
            node.id = :bottom
        end

        #generate top nodes
        xinit = dx / 2
        yinit = dy / 2

        topnodes = [Node([dx * (i-1) + xinit, 
            dy * (j-1) + yinit, 
            z0 + interpolator(dx * (i-1) / xmax, dy * (j-1) /ymax)], 
            :free) for i in 1:nx, j in 1:ny]

        for node in topnodes
            node.id = :top
        end
    else
        #generate nodes for bottom plane
        bottomnodes = [Node([dx * (i-1), dy * (j-1), 0.] .+ base, :free) for i in 1:nx+1, j in 1:ny+1]
        for node in bottomnodes
            node.id = :bottom
        end

        #generate top nodes
        xinit = dx / 2
        yinit = dy / 2

        topnodes = [Node([dx * (i-1) + xinit, 
            dy * (j-1) + yinit, 
            z0 + interpolator(dx * (i-1) / xmax, dy * (j-1) /ymax)], 
            :free) for i in 1:nx, j in 1:ny]

        for node in topnodes
            node.id = :top
        end
    end

    

    #elements
    elements = Vector{AbstractElement{Float64}}()

    #generate bottom horizontal elements
    #parallel to x
    for j = 1:ny+1
        for i = 1:nx
            element = TrussElement(bottomnodes[i,j], bottomnodes[i+1,j], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #parallel to y
    for i = 1:nx+1
        for j = 1:ny
            element = TrussElement(bottomnodes[i,j], bottomnodes[i,j+1], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #generate top horizontal elements
    #parallel to x
    for j = 1:ny
        for i = 1:nx-1
            element = TrussElement(topnodes[i,j], topnodes[i+1,j], section)
            element.id = :top

            push!(elements, element)
        end
    end

    #parallel to y
    for i = 1:nx
        for j = 1:ny-1
            element = TrussElement(topnodes[i,j], topnodes[i,j+1], section)
            element.id = :top

            push!(elements, element)
        end
    end

    #generate web elements
    for i = 1:nx
        for j = 1:ny
            e1 = TrussElement(topnodes[i,j], bottomnodes[i,j], section)
            e2 = TrussElement(topnodes[i,j], bottomnodes[i+1,j], section)
            e3 = TrussElement(topnodes[i,j], bottomnodes[i,j+1], section)
            e4 = TrussElement(topnodes[i,j], bottomnodes[i+1,j+1], section)

            e1.id = e2.id = e3.id = e4.id = :web

            push!(elements, e1, e2, e3, e4)
        end
    end


    #generate node index matrices
    count = 1

    ibottomnodes = zeros(Int64, nx+1, ny+1)
    for j in 1:ny+1
        for i in 1:nx+1
            ibottomnodes[i,j] = count
            count += 1
        end
    end

    itopnodes = zeros(Int64, nx, ny)
    for j in 1:ny
        for i in 1:nx
            itopnodes[i,j] = count
            count += 1
        end
    end

    isquares = [[ibottomnodes[i,j],
        ibottomnodes[i+1,j],
        ibottomnodes[i,j+1],
        ibottomnodes[i+1,j+1]] for i in 1:nx, j in 1:ny]


    ix1 = ibottomnodes[1,:]
    ix2 = ibottomnodes[end,:]
    iy1 = ibottomnodes[:,1]
    iy2 = ibottomnodes[:,end]

    #generate supports
    if support == :corner
        for i in [1, nx]
            for j in [1, ny]
                idset = isquares[i,j]

                for k in idset
                    fixnode!(bottomnodes[k], :pinned)
                    bottomnodes[k].id = :support
                end
            end
        end
    elseif support == :center
        isupport = nx % 2 == 1 ? [Int(ceil(nx / 2))] : Int.([ceil(nx / 2), floor(nx /2)])
        jsupport = ny % 2 == 1 ? [Int(ceil(ny / 2))] : Int.([ceil(ny / 2), floor(ny /2)])

        iset = vcat([vec(isquares[i,j]) for i in isupport, j in jsupport]...)

        for i in iset
            bottomnodes[i].id = :support
            fixnode!(bottomnodes[i], :pinned)
        end
    elseif support == :x
        for i in [ix1; ix2]
            bottomnodes[i].id = :support
            fixnode!(bottomnodes[i], :pinned)
        end

    elseif support == :y
        for i in [iy1; iy2]
            bottomnodes[i].id = :support
            fixnode!(bottomnodes[i], :pinned)
        end

    elseif support == :xy
        for i in [ix1; ix2; iy1; iy2]
            bottomnodes[i].id = :support
            fixnode!(bottomnodes[i], :pinned)
        end
    else
        error("support must be: :corner, :center, :x, :y, :xy")
    end

    flatnodes = [vec(bottomnodes); vec(topnodes)]

    #generate loads
    loads  = [NodeForce(node, load) for node in flatnodes if node.id != :support]

    #assemble
    model = Model(flatnodes, elements, loads)
    solve!(model)

    isupport = findall(model.nodes, :support)




    spaceframe = SpaceFrame(model,
        nx,
        dx,
        ny,
        dy,
        z0,
        section,
        support,
        load,
        base,
        ibottomnodes,
        itopnodes,
        isquares,
        isupport,
        ix1,
        ix2,
        iy1,
        iy2)

    return spaceframe
end