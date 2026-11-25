#=
Multipoint constraints (rigid diaphragms) as a linear map AFTER assembly:
u_free = C·û, K̂ = CᵀKC, F̂ = CᵀF. K̂ is linear in K's entries, so
nzval(K̂) = M·nzval(K) with a sparse M built from C. Element assembly,
force recovery and reactions are untouched (they work on the full u).
=#

"""
    ConstraintMap{T}

The constraint transformation of a processed model with rigid diaphragms
(stored in `AnalysisCache.constraint`).

- `red::Vector{Int}`: free index → reduced index (0 for constrained slave DOFs)
- `C::SparseMatrixCSC{T,Int}`: `n_free × n_red`, `u_free = C·û`
- `csrc`: per stored entry of `C`, `(comp, sign, slave, master)` — `comp = 0`
  is the constant `sign`, else `sign·(X[comp, slave] − X[comp, master])`
  (refreshed from the current node positions at every solve)
- `Khat::SparseMatrixCSC{T,Int}`: reduced stiffness (frozen pattern)
- `mrow, mcol, mq1, mq2`: COO of `M` — entry k adds
  `C.nzval[mq1]·C.nzval[mq2]·K.nzval[mcol]` to `Khat.nzval[mrow]`
- `M::SparseMatrixCSC{T,Int}`: `M` assembled at the current coefficients
  (the pure path's differentiable linear map), `nnz(Khat) × nK`
- `nK::Int`: `nnz(K)`
"""
mutable struct ConstraintMap{T}
    red::Vector{Int}
    C::SparseMatrixCSC{T,Int}
    csrc::Vector{NTuple{4,Int}}
    Khat::SparseMatrixCSC{T,Int}
    mrow::Vector{Int}
    mcol::Vector{Int}
    mq1::Vector{Int}
    mq2::Vector{Int}
    M::SparseMatrixCSC{T,Int}
    nK::Int
end

# global DOFs constrained by a diaphragm: slave dof => [(master dof, comp, sign)],
# coefficient = sign (comp = 0) or sign·(x_slave − x_master)[comp]
function _diaphragm_rules(model::Model)
    rules = Dict{Int,Vector{NTuple{3,Int}}}()
    seen = IdDict{Node,Bool}()
    for d in model.diaphragms
        a = d.axis
        b, c = _inplane_axes(a)
        sm = node_dof_start(d.master)
        for s in d.nodes
            haskey(seen, s) && throw(ArgumentError("node $(s.index) belongs to more than one diaphragm"))
            seen[s] = true
            ss = node_dof_start(s)
            # u_b,s = u_b,m − θ_a,m·r_c ; u_c,s = u_c,m + θ_a,m·r_b ; θ_a,s = θ_a,m
            rules[ss+b] = [(sm + b, 0, 1), (sm + 3 + a, c, -1)]
            rules[ss+c] = [(sm + c, 0, 1), (sm + 3 + a, b, 1)]
            rules[ss+3+a] = [(sm + 3 + a, 0, 1)]
        end
    end
    return rules
end

"""
    build_constraint(model, part, K) -> Union{Nothing, ConstraintMap}

Symbolic construction of the constraint map (once per topology, in
`process!`). Returns `nothing` for models without diaphragms.
"""
function build_constraint(model::Model{T}, part::DofPartition, K::SparseMatrixCSC{T,Int}) where {T}
    isempty(model.diaphragms) && return nothing
    g2f = part.global_to_free
    nfree = length(part.free)
    master_of = Dict{Int,Int}()                    # slave node index → master node index
    for d in model.diaphragms, s in d.nodes
        master_of[s.index] = d.master.index
    end
    rules = _diaphragm_rules(model)
    for (g, _) in rules
        g2f[g] == 0 && g in part.fixed && throw(ArgumentError(
            "a diaphragm slave DOF (global $g) is supported — remove the in-plane support " *
            "or leave the node out of the diaphragm"))
    end
    # reduced numbering: free DOFs that are not constrained slaves
    red = zeros(Int, nfree)
    nred = 0
    for f in 1:nfree
        haskey(rules, part.free[f]) && continue
        nred += 1
        red[f] = nred
    end
    I_ = Int[]; J_ = Int[]; src = NTuple{4,Int}[]
    for f in 1:nfree
        g = part.free[f]
        if haskey(rules, g)
            slave = (g - 1) ÷ 6 + 1
            for (gm, comp, sgn) in rules[g]
                fm = g2f[gm]
                fm == 0 && continue                 # master DOF fixed/inactive: contributes nothing
                push!(I_, f); push!(J_, red[fm]); push!(src, (comp, sgn, slave, master_of[slave]))
            end
        else
            push!(I_, f); push!(J_, red[f]); push!(src, (0, 1, 0, 0))
        end
    end
    C = sparse(I_, J_, ones(T, length(I_)), nfree, nred)
    # map each source entry to its stored position in C (no duplicates by construction)
    pos = _nz_position_lookup(C)
    csrc = Vector{NTuple{4,Int}}(undef, nnz(C))
    for k in eachindex(I_)
        csrc[pos(I_[k], J_[k])] = src[k]
    end
    # rows of C as (col, nz index) lists
    rowent = [Tuple{Int,Int}[] for _ in 1:nfree]
    for j in 1:nred, q in nzrange(C, j)
        push!(rowent[rowvals(C)[q]], (j, q))
    end
    # K̂ pattern and the COO of M
    HI = Int[]; HJ = Int[]
    q1s = Int[]; q2s = Int[]; ps = Int[]
    for jK in 1:nfree, p in nzrange(K, jK)
        iK = rowvals(K)[p]
        for (r, q1) in rowent[iK], (s, q2) in rowent[jK]
            push!(HI, r); push!(HJ, s); push!(ps, p); push!(q1s, q1); push!(q2s, q2)
        end
    end
    Khat = sparse(HI, HJ, zeros(T, length(HI)), nred, nred)
    fill!(nonzeros(Khat), zero(T))
    hpos = _nz_position_lookup(Khat)
    mrow = [hpos(HI[k], HJ[k]) for k in eachindex(HI)]
    cm = ConstraintMap{T}(red, C, csrc, Khat, mrow, ps, q1s, q2s, spzeros(T, nnz(Khat), nnz(K)), nnz(K))
    return cm
end

# refresh C's offset coefficients from node positions (columns of X) and rebuild M
function _refresh_constraint!(cm::ConstraintMap{T}, model::Model{T}) where {T}
    nz = nonzeros(cm.C)
    @inbounds for k in eachindex(nz)
        comp, sgn, s, m = cm.csrc[k]
        nz[k] = comp == 0 ? T(sgn) :
                sgn * (model.nodes[s].position[comp] - model.nodes[m].position[comp])
    end
    vals = similar(nz, length(cm.mrow))
    @inbounds for k in eachindex(vals)
        vals[k] = nz[cm.mq1[k]] * nz[cm.mq2[k]]
    end
    cm.M = sparse(cm.mrow, cm.mcol, vals, nnz(cm.Khat), cm.nK)
    return cm
end

# in-place reduced stiffness values from K's
function _reduce_K!(cm::ConstraintMap{T}, K::SparseMatrixCSC{T,Int}) where {T}
    hz = nonzeros(cm.Khat)
    fill!(hz, zero(T))
    cz = nonzeros(cm.C)
    kz = nonzeros(K)
    @inbounds for k in eachindex(cm.mrow)
        hz[cm.mrow[k]] += cz[cm.mq1[k]] * cz[cm.mq2[k]] * kz[cm.mcol[k]]
    end
    return cm.Khat
end

# pure (differentiable) reduced stiffness: a sparse linear map of K's values
_reduce_K(cm::ConstraintMap, K::SparseMatrixCSC) = sparse_from_pattern(cm.Khat, cm.M * nonzeros(K))
