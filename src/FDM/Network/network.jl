"""
    Network(nodes::Vector{FDMnode}, elements::Vector{FDMelement}, loads::Vector{FDMload})

A FDM network defined by a collection of nodes, elements, and loads.

# Fields
- `nodes::Vector{FDMnode}` collection of nodes in network
- `elements::Vector{FDMelement}` collection of elements in network
- `loads::Vector{FDMload}` collection of external loads in network
- `q::Vector{<:Real}` vector of elemental force densities [n_e × 1]
- `Q::SparseMatrixCSC{Float64, Int64}` diagonalized sparse matrix of `q` [n_e × n_e]
- `C::SparseMatrixCSC{Int64, Int64}` element/node connectivity matrix [n_e × n_n]
- `N::Vector{Int64}` vector of free node indices
- `F::Vector{Int64}` vector of fixed node indices
- `Cn::SparseMatrixCSC{Int64, Int64}` connectivity matrix of free nodes = `C[:, N]`
- `Cf::SparseMatrixCSC{Int64, Int64}` connectivity matrix of fixed nofes = `C[:. F]`
- `P::Matrix{<:Real}` external load matrix [n_n × 3]
- `Pn::Matrix{<:Real}` loads applied to free nodes = `P[N, :]`
- `xyz::Matrix{<:Real}` positions of all nodes [n_n × 3]
"""
mutable struct Network
    nodes::Vector{FDMnode} #vector of nodes
    elements::Vector{FDMelement} #vector of elements
    loads::Vector{FDMload} #vector of loads
    q::Vector{<:Real} #vector of force densities
    Q::SparseMatrixCSC{Float64, Int64} #diagm(q)
    C::SparseMatrixCSC{Int64, Int64} #branch node matrix
    N::Vector{Int64} #fixed node indices
    F::Vector{Int64} #free node indices
    Cn::SparseMatrixCSC{Int64, Int64} #fixed branch node matrix
    Cf::SparseMatrixCSC{Int64, Int64} #free branch node matrix
    P::Matrix{<:Real}#load matrix
    Pn::Matrix{<:Real}#free node load matrix
    xyz::Union{Matrix{Int64}, Matrix{Float64}} #xyz matrix of initial positions
    processed::Bool

    function Network(nodes::Vector{FDMnode}, elements::Vector{FDMelement}, loads::Vector{FDMload})

        network = new(nodes, elements, loads)

        network.processed = false
        return network
    end
end
