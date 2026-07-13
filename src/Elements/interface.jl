"""
    AbstractElement{T<:Real}

Supertype of all structural elements.

An element is **pure definition data** — which nodes it connects, its
section, its end conditions. All analysis state of the legacy library
(cached stiffness, transformation, global DOF ids, forces) lives elsewhere:
the analysis structure in the model's `AnalysisCache`, results in
`LinearResults`.

# The element interface

Every element type implements this contract; the analysis layer is written
against it, so adding an element type (or, later, mass and geometric
stiffness) touches no assembly code.

Topology:
- `nodes(el)` — tuple of connected `Node`s
- `n_internal_dofs(el)` — extra non-nodal DOFs the element requests
  (0 for primitives; used by super-elements like `VariableElement`)
- `ndofs(el)` — total DOF count = `6 × length(nodes(el)) + n_internal_dofs(el)`
- `dof_signature(el)` — `NTuple{ndofs, Bool}`: which of its DOF slots the
  element actually couples stiffness to. Drives global DOF *activity*: a
  slot no element touches never enters the solve. A truss element's
  signature is true only on translations — this is what removes the
  rotational DOFs of truss-only nodes, and the free-torsion singularity of
  hinge-released members, structurally.

Physics (pure kernels of positions and properties — shared by the in-place
and AD assembly paths):
- `stiffness(el, x1, x2)` — element stiffness in GLOBAL coordinates as an
  `SMatrix{n,n}`, `n = ndofs(el)`, with positions passed explicitly

Future hooks (declared now so dynamics/nonlinearity reuse the assembly
machinery; implementations arrive with those analyses):
- `mass(el, x1, x2)`
- `geometric_stiffness(el, x1, x2, N)`
"""
abstract type AbstractElement{T<:Real} end

function nodes end
function dof_signature end
function stiffness end

"""
    n_internal_dofs(el) -> Int

Number of non-nodal (internal) DOFs the element requests from the model's
global DOF space. Zero for primitive elements; super-elements with interior
joints (e.g. `VariableElement`) return `6 × (number of interior joints)`.
"""
n_internal_dofs(::AbstractElement) = 0

"""
    ndofs(el) -> Int

Total number of DOF slots the element maps to: six per connected node plus
any internal DOFs.
"""
ndofs(el::AbstractElement) = 6 * length(nodes(el)) + n_internal_dofs(el)

"""
    mass(el, x1, x2)

Consistent mass matrix hook (global coordinates). Not yet implemented for
any element — arrives with dynamic analysis; declared so the assembly layer
can be written against it today.
"""
mass(el::AbstractElement, x1, x2) =
    error("mass matrix not yet implemented for $(typeof(el)) — dynamics is a post-v1.0 feature")

"""
    geometric_stiffness(el, x1, x2, N)

Geometric (initial-stress) stiffness hook given axial force `N`. Not yet
implemented — arrives with geometric nonlinear analysis.
"""
geometric_stiffness(el::AbstractElement, x1, x2, N) =
    error("geometric stiffness not yet implemented for $(typeof(el)) — nonlinear analysis is a post-v1.0 feature")
