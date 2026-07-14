"""
    CornerDome <: AbstractGenerator

A spaceframe dome supported at four corners

# Fields
- `model::Model{Float64}` structural model
- `nx::Integer` number of bays in x direction
- `dx::Real` x-spacing
- `ny::Integer` number of bays in y direction
- `dy::Real` y-spacing
- `dz::Real` z-offset between top and bottom layers of dome
- `section::AbstractSection` element cross section
"""
struct CornerDome <: AbstractGenerator
    model::Model{Float64}
    nx::Integer
    dx::Real
    ny::Integer
    dz::Real
    section::AbstractSection

    
end

"""
    CornerDome(nx::Integer, dx::Real, ny::Integer, dy::Real, dz::Real, interpolator, section::AbstractSection, load = [0., 0., -10.])

Create a CornerDome.

# Arguments
- `nx::Integer` number of bays in x direction
- `dx::Real` x-spacing
- `ny::Integer` number of bays in y direction
- `dy::Real` y-spacing
- `dz::Real` z-offset between top and bottom layers of dome
- `interpolator` surface height function (see below)
- `section::AbstractSection` element cross section

# `interpolator`
Any callable `interpolator(u, v)` mapping normalized plan coordinates
`u, v ∈ [0, 1]` to the surface Z of the structure — a closure over an
analytic surface, or e.g. an `Interpolations.jl` extrapolation built on a
`(0,1) × (0,1)` grid:

```julia
n = 4
x = range(0,1,n)
y = range(0,1,n)
z = rand(n,n) .* 4500

interpolator = cubic_spline_interpolation((x,y), z)
```
"""
function CornerDome(nx::Integer, dx::Real, ny::Integer, dy::Real, dz::Real, interpolator, section::AbstractSection, load = [0., 0., -10.])

    xmax = dx * nx
    ymax = dy * ny

    #generate nodes for bottom plane
    bottomnodes = [Node([dx * (i-1), 
        dy * (j-1), 
        interpolator(dx * (i-1) / xmax, dy * (j-1) / ymax)], :free, :bottom) for i in 1:nx+1, j in 1:ny+1]
    for node in bottomnodes
        node.id = :bottom
    end

    #generate top nodes
    xinit = dx / 2
    yinit = dy / 2

    topnodes = [Node([dx * (i-1) + xinit, 
        dy * (j-1) + yinit, 
        dz + interpolator(dx * (i-1) / xmax, dy * (j-1) /ymax)], 
        :free) for i in 1:nx, j in 1:ny]

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

    #generate supports
    fixnode!(bottomnodes[1,1], :pinned)
    fixnode!(bottomnodes[1,end], :pinned)
    fixnode!(bottomnodes[end,1], :pinned)
    fixnode!(bottomnodes[end,end], :pinned)

    bottomnodes[1,1].id = bottomnodes[1,end].id = bottomnodes[end,1].id = bottomnodes[end,end].id = :support
    flatnodes = [vec(bottomnodes); vec(topnodes)]

    #generate loads
    loads  = [NodeForce(node, load) for node in flatnodes if node.id != :support]

    #assemble
    truss = Model(flatnodes, elements, loads)
    solve!(truss)

    CornerDome(
        truss,
        nx,
        dx,
        ny,
        dz,
        section
    )

end