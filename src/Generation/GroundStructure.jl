abstract type GroundStructure end

"""
    XGroundStructure(Lx, nx, Ly, ny)

Returns a ground structure topology object.

# Arguments
- `Lx` length in x direction
- `nx` number of bays in x direction
- `Ly` length in y direction
- `ny` number of bays in y direction

# Return
Returns a `XGroundStructure` object.

    struct XGroundStructure <: GroundStructure
        Lx::Float64
        nx::Float64
        dx::Float64
        Ly::Float64
        ny::Float64
        dy::Float64
        xy::Vector{Vector{Float64}}
        Xmatrix::Matrix{Float64}
        Ymatrix::Matrix{Float64}
        igrid::Matrix{Int64}
        ielements::Vector{Vector{Int64}}
        C::SparseMatrixCSC{Int64, Int64}

## Additional Fields
- `igrid` Matrix of node indices with respect to the generated grid
- `Xmatrix` Matrix of nodal x positions with respect to the generated grid
- `Ymatrix` Matrix of nodal y positions with respect to the generated grid
- `ielements` Element start/end nodal indices
- `C` An [nelement × nnode] connectivity matrix
"""
struct XGroundStructure <: GroundStructure
    Lx::Float64
    nx::Float64
    dx::Float64
    Ly::Float64
    ny::Float64
    dy::Float64
    xy::Vector{Vector{Float64}}
    Xmatrix::Matrix{Float64}
    Ymatrix::Matrix{Float64}
    igrid::Matrix{Int64}
    ielements::Vector{Vector{Int64}}
    C::SparseMatrixCSC{Int64, Int64}

    function XGroundStructure(Lx::Float64, nx::Integer, Ly::Float64, ny::Integer)

        x_positions = range(0, Lx, nx)
        y_positions = range(0, Ly, ny)

        dx = Lx / (nx-1)
        dy = Ly / (ny-1)

        xy = Vector{Vector{Float64}}()
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

                push!(xy, [x, y])
                Xmatrix[ix, iy] = x
                Ymatrix[ix, iy] = y

            end
        end

        # make horizontal elements
        ielements = Vector{Vector{Int64}}()

        I = Vector{Int64}()
        J = Vector{Int64}()
        V = Vector{Int64}()

        i_element = 1

        for i = 1:ny
            for j = 1:nx-1
                push!(ielements, [igrid[i,j], igrid[i,j+1]])
                push!(I, i_element, i_element)
                push!(J, igrid[i,j], igrid[i,j+1])
                push!(V, -1, 1)
                i_element += 1
            end
        end

        #vertical elements
        for j = 1:nx
            for i = 1:ny-1
                push!(ielements, [igrid[i,j], igrid[i+1,j]])
                push!(I, i_element, i_element)
                push!(J, igrid[i,j], igrid[i+1,j])
                push!(V, -1, 1)
                i_element += 1
            end
        end

        #diagonal set 1
        for i = 1:ny-1
            for j = 1:nx-1
                push!(ielements, [igrid[i,j], igrid[i+1,j+1]])
                push!(I, i_element, i_element)
                push!(J, igrid[i,j], igrid[i+1,j+1])
                push!(V, -1, 1)
                i_element += 1
            end
        end

        #diagonal set 2
        for i = 2:ny
            for j = 1:nx-1
                push!(ielements, [igrid[i,j], igrid[i-1,j+1]])
                push!(I, i_element, i_element)
                push!(J, igrid[i,j], igrid[i-1,j+1])
                push!(V, -1, 1)
                i_element += 1
            end
        end

        C = sparse(I, J, V)

        new(
            Lx,
            nx,
            dx,
            Ly,
            ny,
            dy,
            xy,
            Xmatrix,
            Ymatrix,
            igrid,
            ielements,
            C
        )
    end
end

"""
    DenseGroundStructure(Lx, nx, Ly, ny)

Returns a ground structure topology object.

# Arguments
- `Lx` length in x direction
- `nx` number of bays in x direction
- `Ly` length in y direction
- `ny` number of bays in y direction

# Return
Returns a `DenseGroundStructure` object.

    struct DenseGroundStructure <: GroundStructure
        Lx::Float64
        nx::Float64
        dx::Float64
        Ly::Float64
        ny::Float64
        dy::Float64
        xy::Vector{Vector{Float64}}
        Xmatrix::Matrix{Float64}
        Ymatrix::Matrix{Float64}
        igrid::Matrix{Int64}
        ielements::Vector{Vector{Int64}}
        C::SparseMatrixCSC{Int64, Int64}

## Additional Fields
- `igrid` Matrix of node indices with respect to the generated grid
- `Xmatrix` Matrix of nodal x positions with respect to the generated grid
- `Ymatrix` Matrix of nodal y positions with respect to the generated grid
- `ielements` Element start/end nodal indices
- `C` An [nelement × nnode] connectivity matrix
"""
struct DenseGroundStructure <: GroundStructure
    Lx::Float64
    nx::Float64
    dx::Float64
    Ly::Float64
    ny::Float64
    dy::Float64
    xy::Vector{Vector{Float64}}
    Xmatrix::Matrix{Float64}
    Ymatrix::Matrix{Float64}
    igrid::Matrix{Int64}
    ielements::Vector{Vector{Int64}}
    C::SparseMatrixCSC{Int64, Int64}

    function DenseGroundStructure(Lx::Float64, nx::Integer, Ly::Float64, ny::Integer)

        x_positions = range(0, Lx, nx)
        y_positions = range(0, Ly, ny)

        dx = Lx / (nx-1)
        dy = Ly / (ny-1)

        xy = Vector{Vector{Float64}}()
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

                push!(xy, [x, y])
                Xmatrix[ix, iy] = x
                Ymatrix[ix, iy] = y

            end
        end

        # make horizontal elements
        ielements = Vector{Vector{Int64}}()

        index_vector = igrid[:]
        n_nodes = length(index_vector)

        I = Vector{Int64}()
        J = Vector{Int64}()
        V = Vector{Int64}()

        i_element = 1
        for i = 1:n_nodes-1
            for j = i+1:n_nodes
                push!(ielements, [i,j])
                
                push!(I, i_element, i_element)
                push!(J, i, j)
                push!(V, -1, 1)

                i_element += 1
            end
        end

        C = sparse(I, J, V)

        new(
            Lx,
            nx,
            dx,
            Ly,
            ny,
            dy,
            xy,
            Xmatrix,
            Ymatrix,
            igrid,
            ielements,
            C
        )

    end
end

"""
    BoundedGroundStructure(Lx, nx, Ly, ny, ub)

Returns a ground structure topology object.

# Arguments
- `Lx` length in x direction
- `nx` number of bays in x direction
- `Ly` length in y direction
- `ny` number of bays in y direction
- `ub` node-node length threshold to generate a connection

# Return
Returns a `BoundedGroundStructure` object.

    struct BoundedGroundStructure <: GroundStructure
        Lx::Float64
        nx::Float64
        dx::Float64
        Ly::Float64
        ny::Float64
        dy::Float64
        xy::Vector{Vector{Float64}}
        Xmatrix::Matrix{Float64}
        Ymatrix::Matrix{Float64}
        igrid::Matrix{Int64}
        ielements::Vector{Vector{Int64}}
        C::SparseMatrixCSC{Int64, Int64}

## Additional Fields
- `igrid` Matrix of node indices with respect to the generated grid
- `Xmatrix` Matrix of nodal x positions with respect to the generated grid
- `Ymatrix` Matrix of nodal y positions with respect to the generated grid
- `ielements` Element start/end nodal indices
- `C` An [nelement × nnode] connectivity matrix
"""
struct BoundedGroundStructure <: GroundStructure
    Lx::Float64
    nx::Float64
    dx::Float64
    Ly::Float64
    ny::Float64
    dy::Float64
    ub::Float64
    xy::Vector{Vector{Float64}}
    Xmatrix::Matrix{Float64}
    Ymatrix::Matrix{Float64}
    igrid::Matrix{Int64}
    ielements::Vector{Vector{Int64}}
    C::SparseMatrixCSC{Int64, Int64}

    function BoundedGroundStructure(Lx::Float64, nx::Integer, Ly::Float64, ny::Integer, ub::Float64)

        x_positions = range(0, Lx, nx)
        y_positions = range(0, Ly, ny)

        dx = Lx / (nx-1)
        dy = Ly / (ny-1)

        @assert dx ≤ ub && dy ≤ ub "Upper bound is smaller than grid size, reduce bound or increase density"

        xy = Vector{Vector{Float64}}()
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

                push!(xy, [x, y])
                Xmatrix[ix, iy] = x
                Ymatrix[ix, iy] = y

            end
        end

        # make horizontal elements
        i_padded_col = Int(floor(ub / dx))
        i_padded_row = Int(floor(ub / dy))

        padded_index_matrix = [
            zeros(Int64, i_padded_row, nx + 2i_padded_col);
            zeros(Int64, ny, i_padded_col) igrid zeros(Int64, ny, i_padded_col);
            zeros(Int64, i_padded_row, nx + 2i_padded_col)
            ]


        ielements = Vector{Vector{Int64}}()
        I = Vector{Int64}()
        J = Vector{Int64}()
        V = Vector{Int64}()

        i_element = 1
        
        for i = 1+i_padded_row:ny+i_padded_row
            for j = 1+i_padded_col:nx+i_padded_col

                i_node = padded_index_matrix[i,j]

                threshold_indices = padded_index_matrix[i-i_padded_row:i+i_padded_row, j-i_padded_col:j+i_padded_col]

                i_valid = [index for index in threshold_indices if index > 0 && index != i_node]

                distances = [norm(xy[index] - xy[i_node]) for index in i_valid]

                i_to_connect = i_valid[findall(distances .≤ ub)]

                for index in i_to_connect
                    push!(ielements, [i_node, index])
                    push!(I, i_element, i_element)
                    push!(J, i_node, index)
                    push!(V, -1, 1)
                    i_element += 1
                end

                padded_index_matrix[i,j] = 0

            end
        end

        C = sparse(I, J, V)
        
        new(
            Lx,
            nx,
            dx,
            Ly,
            ny,
            dy,
            ub,
            xy,
            Xmatrix,
            Ymatrix,
            igrid,
            ielements,
            C
        )

    end
end