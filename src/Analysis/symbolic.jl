"""
    ElementGroup{E<:AbstractElement}

Elements of one concrete type, gathered so numeric assembly runs
type-stably behind a function barrier (the model's element vector is
heterogeneous by design).

Per element (parallel vectors, indexed within the group):
- `elements::Vector{E}`
- `model_indices`: position of each element in `model.elements`
- `slots`: the element's *active* local DOF slots (from its
  [`dof_signature`](@ref)) — only these are assembled
- `gdofs`: the global DOF index of each active slot
- `nzmap`: for each active `(i, j)` slot pair (column-major over
  `length(slots)²`), the linear index into the frozen `K.nzval` where that
  stiffness entry accumulates, or `0` if the row or column DOF is not free
"""
struct ElementGroup{E<:AbstractElement}
    elements::Vector{E}
    model_indices::Vector{Int}
    slots::Vector{Vector{Int}}
    gdofs::Vector{Vector{Int}}
    nzmap::Vector{Vector{Int}}
end

"""
    AnalysisCache{T}

The analysis structure derived from a `Model` by [`process!`](@ref) —
everything the solver needs that depends only on *topology*, built once and
reused across solves:

- `partition::DofPartition`: free/fixed/inactive classification
- `groups::Vector{ElementGroup}`: type-grouped elements with frozen scatter
  maps
- `K::SparseMatrixCSC{T,Int}`: the free×free stiffness matrix with a FROZEN
  sparsity pattern — numeric assembly only rewrites `nonzeros(K)`; no
  full-space matrix is ever formed
- `spring_nz::Vector{Tuple{Int,T}}`: (nzval position, stiffness) pairs for
  nodal springs on free DOFs
- `P`, `Pf::Vector{T}`: full-space nodal-load and fixed-end-force vectors
- `q_local::Vector{Vector{T}}`: per-element accumulated LOCAL condensed
  fixed-end forces (needed for element force recovery)
- `factorization`: cached factorization of `K` (set on first solve; numeric
  refactorization reuses the symbolic analysis since the pattern is frozen)

Geometry/section/load *values* may change freely between solves; adding or
removing nodes, elements, springs, or changing end conditions in a way that
alters DOF activity requires re-processing (`process!`).
"""
mutable struct AnalysisCache{T}
    partition::DofPartition
    groups::Vector{ElementGroup}
    K::SparseMatrixCSC{T,Int}
    spring_nz::Vector{Tuple{Int,T}}
    P::Vector{T}
    Pf::Vector{T}
    q_local::Vector{Vector{T}}
    factorization::Any
end

# active slots of an element (positions 1:12 where the signature is true)
_active_slots(el::AbstractElement) = findall(collect(dof_signature(el)))

"""
    build_cache(model) -> AnalysisCache

The symbolic assembly pass: partition DOFs, group elements by concrete
type, build the free×free sparsity pattern from all active-slot pairs (plus
spring diagonals), freeze it, and invert it into per-element `nzmap`s so
numeric assembly is pure scatter with zero allocation.
"""
function build_cache(model::Model{T}) where {T}
    part = partition_dofs(model)
    g2f = part.global_to_free
    nfree = length(part.free)

    # group elements by concrete type
    bytype = Dict{DataType,Vector{Int}}()
    for (i, el) in enumerate(model.elements)
        push!(get!(bytype, typeof(el), Int[]), i)
    end

    # collect the pattern: every free×free active-slot pair of every element
    I_ = Int[]
    J_ = Int[]
    protogroups = Vector{Tuple{DataType,Vector{Int},Vector{Vector{Int}},Vector{Vector{Int}}}}()
    for (E, idxs) in sort!(collect(bytype); by=p -> string(p.first))
        slots_all = Vector{Vector{Int}}()
        gdofs_all = Vector{Vector{Int}}()
        for i in idxs
            el = model.elements[i]
            slots = _active_slots(el)
            g = element_global_dofs(el)[slots]
            push!(slots_all, slots)
            push!(gdofs_all, g)
            for a in eachindex(g), b in eachindex(g)
                fa = g2f[g[a]]
                fb = g2f[g[b]]
                if fa > 0 && fb > 0
                    push!(I_, fa)
                    push!(J_, fb)
                end
            end
        end
        push!(protogroups, (E, idxs, slots_all, gdofs_all))
    end
    for sp in model.springs
        s = node_dof_start(sp.node)
        for i in 1:6
            f = g2f[s+i]
            if sp.stiffness[i] > 0 && f > 0
                push!(I_, f)
                push!(J_, f)
            end
        end
    end

    K = sparse(I_, J_, zeros(T, length(I_)), nfree, nfree)
    fill!(nonzeros(K), zero(T))

    # invert the pattern: per-element nzval positions (the legacy AsapOptim
    # `all_inz` hack, now a first-class product of symbolic assembly)
    nzpos = _nz_position_lookup(K)
    groups = ElementGroup[]
    for (E, idxs, slots_all, gdofs_all) in protogroups
        nzmaps = Vector{Vector{Int}}()
        for g in gdofs_all
            n = length(g)
            m = Vector{Int}(undef, n * n)
            p = 0
            for b in 1:n, a in 1:n           # column-major over (a, b)
                p += 1
                fa = g2f[g[a]]
                fb = g2f[g[b]]
                m[p] = (fa > 0 && fb > 0) ? nzpos(fa, fb) : 0
            end
            push!(nzmaps, m)
        end
        els = Vector{E}([model.elements[i] for i in idxs])
        push!(groups, ElementGroup{E}(els, idxs, slots_all, gdofs_all, nzmaps))
    end

    spring_nz = Tuple{Int,T}[]
    for sp in model.springs
        s = node_dof_start(sp.node)
        for i in 1:6
            f = g2f[s+i]
            if sp.stiffness[i] > 0 && f > 0
                push!(spring_nz, (nzpos(f, f), T(sp.stiffness[i])))
            end
        end
    end

    return AnalysisCache{T}(part, groups, K, spring_nz,
        zeros(T, part.n_global), zeros(T, part.n_global),
        [zeros(T, 12) for _ in model.elements], nothing)
end

# closure looking up the nzval position of entry (i, j) in a CSC matrix
function _nz_position_lookup(K::SparseMatrixCSC)
    colptr = K.colptr
    rowval = K.rowval
    return function (i::Int, j::Int)
        r = searchsorted(view(rowval, colptr[j]:colptr[j+1]-1), i)
        isempty(r) && error("entry ($i, $j) not in the frozen sparsity pattern")
        return colptr[j] + first(r) - 1
    end
end
