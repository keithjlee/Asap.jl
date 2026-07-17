#=
The FDM analysis structure — the same three-layer separation as the frame
core: Network (definition) / NetworkCache (per-topology structure) /
NetworkResults (per-solve results).
=#

"""
    AxisSystem{T}

One force-density system `D · x = b`: the free/fixed node partition of one
or more coordinate axes (axes with identical partitions share a system, so
uniform networks hold a single system solving all three coordinates as one
multi-RHS problem).

`D = CₙᵀQCₙ` is the free-node force-density Laplacian: its sparsity pattern
depends only on topology, so it is FROZEN here and numeric assembly is a
pure scatter of the current force densities through `nzmap` —
`D[i,i] = Σq` at node i, `D[i,j] = −Σq` of members connecting i and j.

# Fields
- `axes::Vector{Int}`: which coordinates (1 = x, 2 = y, 3 = z) this system solves
- `Na`, `Fa::Vector{Int}`: free/fixed node indices of this partition
- `D::SparseMatrixCSC{T,Int}`: the frozen-pattern system matrix
- `nzmap::Vector{NTuple{4,Int}}`: per ELEMENT, the `nzval` positions of its
  (ii, jj, ij, ji) contributions (0 = that end is fixed here)
- `bnd::Vector{Tuple{Int,Int,Int}}`: (element index, local row, fixed node
  index) triples — boundary members contribute `q · x_fixed` to the RHS
- `fact`: the seam [`FactorizationCache`](@ref) — reused across solves
  (numeric-only refactorization; any solver backend)
"""
mutable struct AxisSystem{T}
    axes::Vector{Int}
    Na::Vector{Int}
    Fa::Vector{Int}
    D::SparseMatrixCSC{T,Int}
    nzmap::Vector{NTuple{4,Int}}
    bnd::Vector{Tuple{Int,Int,Int}}
    fact::Any
end

"""
    NetworkCache{T}

Everything the FDM solver derives from TOPOLOGY, built once by
[`process!`](@ref) and reused across solves: connectivity, node partitions,
and one frozen-pattern [`AxisSystem`](@ref) per distinct axis partition.

Force densities, node positions (anchors included), and loads are read
FRESH from the network at every `solve!` — change them freely; only
adding/removing nodes/elements/loads or changing fixity requires
re-processing.

# Fields
- `C::SparseMatrixCSC{Int,Int}`: element/node incidence (−1 start, +1 end)
- `N`, `F::Vector{Int}`: fully-free / remaining node indices (uniform split)
- `Naxis`, `Faxis::Vector{Vector{Int}}`: per-axis partitions
- `mixed::Bool`: any node with non-uniform fixity?
- `systems::Vector{AxisSystem{T}}`: the solvable systems
- `xyz::Matrix{T}`: [n × 3] working geometry (refreshed from nodes each solve)
- `P::Matrix{T}`: [n × 3] load matrix (refreshed from loads each solve)
"""
mutable struct NetworkCache{T}
    C::SparseMatrixCSC{Int,Int}
    N::Vector{Int}
    F::Vector{Int}
    Naxis::Vector{Vector{Int}}
    Faxis::Vector{Vector{Int}}
    mixed::Bool
    systems::Vector{AxisSystem{T}}
    xyz::Matrix{T}
    P::Matrix{T}
end

"""
    process!(network::Network) -> network

Build the network's analysis structure: assign indices, partition nodes per
axis, and freeze one force-density pattern per distinct partition. Call
once per topology — geometry, force-density, and load VALUES are picked up
fresh by every `solve!`.
"""
function process!(network::Network{T}) where {T}
    nodes = network.nodes
    elements = network.elements
    nn = Base.length(nodes)
    ne = Base.length(elements)

    for (i, node) in enumerate(nodes)
        node.index = i
    end
    for (i, el) in enumerate(elements)
        el.index = i
        el.iStart = el.pStart.index
        el.iEnd = el.pEnd.index
    end

    #incidence
    C = spzeros(Int, ne, nn)
    for (i, el) in enumerate(elements)
        C[i, el.iStart] = -1
        C[i, el.iEnd] = 1
    end

    #partitions
    N = findall(n -> all(n.fixity), nodes)
    F = findall(n -> !all(n.fixity), nodes)
    Naxis = [findall(n -> n.fixity[a], nodes) for a in 1:3]
    Faxis = [findall(n -> !n.fixity[a], nodes) for a in 1:3]
    mixed = any(n -> any(n.fixity) != all(n.fixity), nodes)

    #group axes with identical partitions into shared systems
    systems = AxisSystem{T}[]
    claimed = falses(3)
    for a in 1:3
        claimed[a] && continue
        axes = [b for b in a:3 if Naxis[b] == Naxis[a]]
        claimed[axes] .= true
        isempty(Naxis[a]) && continue
        push!(systems, _axis_system(T, elements, Naxis[a], Faxis[a], axes, nn))
    end

    network.cache = NetworkCache{T}(C, N, F, Naxis, Faxis, mixed, systems,
        zeros(T, nn, 3), zeros(T, nn, 3))
    network.results = nothing
    return network
end

#build one frozen-pattern force-density system for a node partition
function _axis_system(::Type{T}, elements, Na::Vector{Int}, Fa::Vector{Int},
    axes::Vector{Int}, nn::Int) where {T}
    g2f = zeros(Int, nn)
    for (f, g) in enumerate(Na)
        g2f[g] = f
    end

    I_ = Int[]
    J_ = Int[]
    for el in elements
        fi = g2f[el.iStart]
        fj = g2f[el.iEnd]
        if fi > 0
            push!(I_, fi); push!(J_, fi)
        end
        if fj > 0
            push!(I_, fj); push!(J_, fj)
        end
        if fi > 0 && fj > 0
            push!(I_, fi); push!(J_, fj)
            push!(I_, fj); push!(J_, fi)
        end
    end
    nfree = Base.length(Na)
    D = sparse(I_, J_, zeros(T, Base.length(I_)), nfree, nfree)

    nzpos = _nz_position_lookup(D)
    nzmap = Vector{NTuple{4,Int}}(undef, Base.length(elements))
    bnd = Tuple{Int,Int,Int}[]
    for (e, el) in enumerate(elements)
        fi = g2f[el.iStart]
        fj = g2f[el.iEnd]
        nzmap[e] = (fi > 0 ? nzpos(fi, fi) : 0,
            fj > 0 ? nzpos(fj, fj) : 0,
            fi > 0 && fj > 0 ? nzpos(fi, fj) : 0,
            fi > 0 && fj > 0 ? nzpos(fj, fi) : 0)
        fi > 0 && fj == 0 && push!(bnd, (e, fi, el.iEnd))
        fj > 0 && fi == 0 && push!(bnd, (e, fj, el.iStart))
    end

    return AxisSystem{T}(axes, Na, Fa, D, nzmap, bnd, nothing)
end
