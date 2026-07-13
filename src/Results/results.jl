"""
    LinearResults{T}

Results of a linear static solve — displacements, reactions, element end
forces, and compliance. Results live HERE, not on nodes/elements: this
keeps definition objects pure and lets the result scalar type differ from
the model's (e.g. dual numbers flowing through the AD path).

Access through the functions below rather than raw fields where possible —
they handle global-DOF bookkeeping for you.

# Fields
- `u::Vector{T}`: global displacement vector, full DOF space (fixed and
  inactive slots are zero) [length, rad]
- `reactions::Vector{T}`: support reactions, full DOF space (nonzero only
  at fixed DOFs) [force, force·length]
- `element_forces::Vector{Vector{T}}`: per element, the LOCAL end-force
  12-vector `[N₁, Vy₁, Vz₁, T₁, My₁, Mz₁, N₂, Vy₂, Vz₂, T₂, My₂, Mz₂]` —
  forces the element exerts on its ends, in element local coordinates
- `compliance::T`: external work `uᵀF` — the standard stiffness objective

# Accessors
- [`displacement`](@ref)`(res, node) -> SVector{6}`
- [`reaction`](@ref)`(res, node) -> SVector{6}`
- [`element_forces`](@ref)`(res, el) -> Vector{T}`
- [`axial_force`](@ref)`(res, el) -> T`
"""
struct LinearResults{T}
    u::Vector{T}
    reactions::Vector{T}
    element_forces::Vector{Vector{T}}
    compliance::T
end

"""
    displacement(res::LinearResults, node) -> SVector{6}

The node's displacements `(ux, uy, uz, θx, θy, θz)` in global coordinates.
"""
function displacement(res::LinearResults{T}, node::Node) where {T}
    s = node_dof_start(node)
    return SVector{6,T}(res.u[s+1], res.u[s+2], res.u[s+3], res.u[s+4], res.u[s+5], res.u[s+6])
end

"""
    reaction(res::LinearResults, node) -> SVector{6}

Support reactions `(Fx, Fy, Fz, Mx, My, Mz)` at the node in global
coordinates. Zero at unsupported nodes.
"""
function reaction(res::LinearResults{T}, node::Node) where {T}
    s = node_dof_start(node)
    return SVector{6,T}(res.reactions[s+1], res.reactions[s+2], res.reactions[s+3],
        res.reactions[s+4], res.reactions[s+5], res.reactions[s+6])
end

"""
    element_forces(res::LinearResults, el) -> Vector{T}

The element's local end-force vector (see [`LinearResults`](@ref) for the
component ordering). Requires the element's `index` (assigned by
`process!`).
"""
element_forces(res::LinearResults, el::AbstractElement) = res.element_forces[el.index]

"""
    axial_force(res::LinearResults, el) -> T

The element's axial force, tension-positive [force]. Uniform across element
types (slot 7 of the local end-force vector — the axial action at the end
node, which equals the member force for a two-node element).
"""
axial_force(res::LinearResults, el::AbstractElement) = element_forces(res, el)[7]
