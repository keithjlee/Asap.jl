"""
    GridFrame(Lx, nx, Ly, ny, section; load = [0., 0., -1.], support = :corner, support_type = :pinned)

Generate a grid-topology Asap model.

# Arguments
- `Lx::Real` total length in x direction
- `nx::Integer` number of bays in x direction
- `Ly::Real` total length in y direction
- `ny::Integer` number of bays in y direction
- `section::Section` cross section assigned to elements

# Optional arguments
- `load = [0, 0, 1]` force applied to all free nodes
- `support = :corner` nodes to fix. Choose from: :corner, :x, :y, :xy
- `support_type = :pinned` boundary condition on support nodes. See fixDict for valud entries

# Return
Returns a `GridFrame` structure.

    struct GridFrame <: AbstractGenerator
        model::Model
        nx::Integer
        dx::Real
        ny::Integer
        dy::Real
        igrid::Matrix{Int64}

## Fields
Along with bay spacing fields, additional fields are included in the returned `GridNetwork` struct:
- `network::Network` FDM network constructed from parameters
- `igrid::Matrix{Int64}` Matrix of nodal indices in relation to the generated grid


"""
struct GridFrame <: AbstractGenerator
    model::Model{Float64}
    nx::Integer
    dx::Real
    ny::Integer
    dy::Real
    igrid::Matrix{Int64}

    function GridFrame(Lx::Real, nx::Integer, Ly::Real, ny::Integer, section::Section; load = [0., 0., -1.], support = :corner, support_type = :pinned)

        @assert in(support, [:corner, :x, :y, :xy])

        x_positions = range(0, Lx, nx)
        y_positions = range(0, Ly, ny)

        dx = Lx / (nx-1)
        dy = Ly / (ny-1)

        xyz = Vector{Vector{Float64}}()
        Xmatrix = zeros(ny, nx)
        Ymatrix = zeros(ny, nx)
        igrid = zeros(Int64, ny, nx)

        index = 1
        for iy = 1:nx
            for ix = 1:ny

                x = x_positions[iy]
                y = y_positions[ix]

                igrid[ix, iy] = index
                index += 1

                push!(xyz, [x, y, 0.])
                Xmatrix[ix, iy] = x
                Ymatrix[ix, iy] = y

            end
        end

        if support == :corner
            support_indices = [igrid[1, 1], igrid[ny, 1], igrid[1, nx], igrid[ny, nx]]
        elseif support == :x
            support_indices = igrid[[1, ny], :]
        elseif support == :y
            support_indices = igrid[:, [1, nx]]
        else
            support_indices = [igrid[[1, ny], :][:]; igrid[2:ny-1, [1, nx]][:]]
        end

        #make nodes
        nodes = [Node(pos, :free, :free) for pos in xyz]

        #make support nodes
        for node in nodes[support_indices]
            fixnode!(node, support_type)
            node.id = :support
        end

        #make elements
        elements = Vector{AbstractElement{Float64}}()

        #horizontal elements
        for i = 1:ny
            for j = 1:nx-1
                index = [igrid[i,j], igrid[i,j+1]]
                push!(elements, FrameElement(nodes[index]..., section))
            end
        end

        #vertical elements
        for j = 1:nx
            for i = 1:ny-1
                index = [igrid[i,j], igrid[i+1,j]]
                push!(elements, FrameElement(nodes[index]..., section)) 
            end
        end

        #loads
        loads = [NodeForce(node, load) for node in nodes[:free]]

        #assemble
        model = Model(nodes, elements, loads)
        solve!(model)

        new(
            model,
            nx,
            dx,
            ny,
            dy,
            igrid
        )
    end

end