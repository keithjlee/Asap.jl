"""
    solve!(network::Network; reprocess = false, solver = nothing)

Find the network's equilibrium geometry and update the node positions in
place (free axes only — fixed coordinates are boundary data). Forces and
reactions land in `network.results`.

Every solve reads the CURRENT force densities, anchor positions, and loads;
the force-density pattern is frozen per topology, so assembly is a pure
scatter and the factorization reuses its symbolic analysis (numeric-only
refactorization — same solver seam as the frame core, so `solver` accepts
any LinearSolve algorithm when LinearSolve.jl is loaded).

FDM equilibrium is separable per coordinate: axes sharing a free/fixed
partition solve as one multi-RHS system (uniform networks: a single system
for x, y, z); per-axis fixity solves each partition independently.
"""
function solve!(network::Network{T}; reprocess::Bool=false, solver=nothing) where {T}
    (network.cache === nothing || reprocess) && process!(network)
    cache = network.cache::NetworkCache{T}

    _refresh_state!(cache, network)

    for sys in cache.systems
        _solve_system!(cache, sys, network, solver)
    end

    #write equilibrium coordinates back onto the nodes (free axes only)
    for (i, node) in enumerate(network.nodes)
        any(node.fixity) || continue
        node.position = SVector{3,T}(ntuple(a ->
            node.fixity[a] ? cache.xyz[i, a] : node.position[a], 3))
    end

    network.results = _network_results(network)
    return network
end

#pull the CURRENT geometry and loads into the cache work arrays
function _refresh_state!(cache::NetworkCache{T}, network::Network{T}) where {T}
    for (i, node) in enumerate(network.nodes)
        cache.xyz[i, :] = node.position
    end
    fill!(cache.P, zero(T))
    for load in network.loads
        cache.P[load.point.index, :] += load.force
    end
    return cache
end

#assemble D (pure scatter of current q) + RHS, refactorize, backsolve
function _solve_system!(cache::NetworkCache{T}, sys::AxisSystem{T},
    network::Network{T}, solver) where {T}
    els = network.elements
    nz = nonzeros(sys.D)
    fill!(nz, zero(T))
    @inbounds for (e, m) in enumerate(sys.nzmap)
        q = els[e].q
        m[1] > 0 && (nz[m[1]] += q)
        m[2] > 0 && (nz[m[2]] += q)
        m[3] > 0 && (nz[m[3]] -= q)
        m[4] > 0 && (nz[m[4]] -= q)
    end

    B = cache.P[sys.Na, sys.axes]
    @inbounds for (e, row, jg) in sys.bnd      # boundary members: q·x_anchor
        q = els[e].q
        for (c, a) in enumerate(sys.axes)
            B[row, c] += q * cache.xyz[jg, a]
        end
    end

    fc = sys.fact
    if fc isa FactorizationCache && (solver === nothing || solver === fc.solver)
        _refactorize!(fc.solver, fc, sys.D)
    else
        fc = _factorize(solver, sys.D)
        sys.fact = fc
    end
    X = _backsolve(fc.solver, fc, B)

    cache.xyz[sys.Na, sys.axes] = X
    return sys
end

"""
    solve(network::Network, q::Vector{<:Real}) -> xyz

Equilibrium positions for a NEW vector of force densities without mutating
the network (same per-axis semantics as [`solve!`](@ref); fresh
factorization, default backend).
"""
function solve(network::Network{T}, q::Vector{<:Real}) where {T}
    @assert Base.length(q) == Base.length(network.elements) "q and elements must have equal length"
    network.cache === nothing && process!(network)
    cache = network.cache::NetworkCache{T}
    _refresh_state!(cache, network)

    xyzout = copy(cache.xyz)
    for sys in cache.systems
        D = copy(sys.D)
        nz = nonzeros(D)
        fill!(nz, zero(T))
        @inbounds for (e, m) in enumerate(sys.nzmap)
            m[1] > 0 && (nz[m[1]] += q[e])
            m[2] > 0 && (nz[m[2]] += q[e])
            m[3] > 0 && (nz[m[3]] -= q[e])
            m[4] > 0 && (nz[m[4]] -= q[e])
        end
        B = cache.P[sys.Na, sys.axes]
        @inbounds for (e, row, jg) in sys.bnd
            for (c, a) in enumerate(sys.axes)
                B[row, c] += q[e] * cache.xyz[jg, a]
            end
        end
        xyzout[sys.Na, sys.axes] = _factorize(nothing, D) \ B
    end
    return xyzout
end
