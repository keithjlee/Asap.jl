[![DOI](https://zenodo.org/badge/426740094.svg)](https://zenodo.org/doi/10.5281/zenodo.10581559)

![](READMEassets/forces-axo.png)

# Asap.jl
Asap is...
- the anti-SAP (2000)
- results as Soon As Possible
- another Structural Analysis Package

Designed first-and-foremost for information-rich data structures and ease of querying, but always with performance in mind.

See also: [AsapToolkit](https://github.com/keithjlee/AsapToolkit), [AsapOptim](https://github.com/keithjlee/AsapOptim), [AsapHarmonics](https://github.com/keithjlee/AsapHarmonics).

# Installation
Asap.jl is now a registered Julia package. Install through package mode in the REPL:
```julia
pkg> add Asap
```
or
```julia
using Pkg
Pkg.Add("Asap")
```

## Citing
When using or extending this software for research purposes, please cite using the following:

### Bibtex
```
@software{lee_2024_10581560,
  author       = {Lee, Keith Janghyun},
  title        = {Asap.jl},
  month        = jan,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v0.1},
  doi          = {10.5281/zenodo.10581560},
  url          = {https://doi.org/10.5281/zenodo.10581560}
}
```

### Other styles
Or find a pre-written citation in the style of your choice [here](https://zenodo.org/records/10724610) (see the Citation box on the right side). E.g., for APA:
```
Lee, K. J. (2024). Asap.jl (v0.1). Zenodo. https://doi.org/10.5281/zenodo.10581560
```

# Extensions, Related packages
See [AsapToolkit.jl](https://github.com/keithjlee/AsapToolkit) for even more utility and post-processing functions.

# Usage
A structural model is defined by:
```julia
model = Model(nodes, elements, loads)
```
and solved via:
```julia
solve!(model)
```
Which finds the unknown nodal displacement field, $u = S^{-1}(P-P_f)$ where:
- $S$ is the global stiffness matrix (often called $K$)
- $P$ is the global external load vector
- $P_f$ is the fixed end forces induced by element loads (such as a line load on an element)

## `Node`
```julia
mutable struct Node <: AbstractNode
    position::Vector{Float64}
    dof::Vector{Bool}
    nodeID::Int64
    globalID::Vector{Int64}
    reaction::Vector{Float64}
    displacement::Vector{Float64}
    id::Symbol
end
```
We begin with the primary information carrier for structural analysis: nodes with *n* independent degrees of freedom (DOF). They are defined by a spatial position in $\mathbb{R}^3$ as well as a vector of booleans that indicate which DOFs are free to move under load, in order: $T_x, T_y, T_z, R_x, R_y, R_z$ where $T$ is a translational DOF and $R$ is a rotational DOF. E.g.:
```julia
node = Node([0, 15.5, 12.0], [false, false, false, true, true, true])
```
This defines a node at $x = 0; y = 15.5; z = 12$ with a *pinned* support (i.e., translational DOFs are fixed, but rotational DOFs are not).

Some common boundary conditions are provided to you as symbols to use in the constructor:
```julia
node = Node([0.,0.,0.], :free) # all free DOFs
node = Node([0.,0.,0.], :fixed) # all fixed DOFs
node = Node([0.,0.,0.], :xfree) # all DOFs are fixed except Tx. You can also do :yfree or :zfree
node = Node([0.,0.,0.], :xfixed) # all DOFs are free except Tx. You can also do :yfixed or :zfixed
```

Nodes can also include an optional identifier represented as a symbol:
```julia
pin_support = Node(zeros(3), :pinned, :pinsupport)
roller_support = Node([15.1, 0, 0], :xfree, :rollersupport)
free_nodes = [Node(rand(3), :free, :freenodes) for _ = 1:10]

nodes = [pin_support; roller_support; free_nodes]
```

This allows you to index into a vector of nodes using the identifier:
```julia
all_free_nodes = nodes[:freenodes] #returns a vector of nodes with the :freenode identifier
```

Or find the indices of nodes in a vector of nodes that have a given id:
```
i_free_nodes = findall(nodes, :freenodes)
```

## `Element`
```julia
mutable struct Element{R<:Release} <: FrameElement{R}
    section::Section #cross section
    nodeStart::Node #start node
    nodeEnd::Node #end position
    elementID::Int64
    globalID::Vector{Int64} #element global DOFs
    length::Float64 #length of element
    K::Matrix{Float64} # stiffness matrix in GCS
    Q::Vector{Float64} # fixed end forces in GCS
    R::Matrix{Float64} # transformation matrix
    Ψ::Float64 #roll angle
    LCS::Vector{Vector{Float64}} #local coordinate frame (X, y, z)
    forces::Vector{Float64} #elemental forces in LCS
    id::Symbol #optional identifier
end
```
Elements are defined by their start and end nodes, a cross-section, and an optional identifier.

### `Section`
A section defines the mechanical and material properties of an element:
```julia
ibeam_section = Section(A, E, G, Ix, Iy, J)
```
where:
- `A`: area
- `E`: Young's Modulus
- `G`: Shear Modulus
- `Ix`: Moment of inertia in strong axis (often denoted as `Iz` in other FEA programs)
- `Iy`: Moment of inertia in weak axis
- `J`: Torsional constant

**It is up to you to ensure unit consistency**.

You can then define an element via:
```julia
element = Element(pin_support, roller_support, ibeam_section)
element_with_id = Element(pin_support, rand(free_nodes), ibeam_section, :randomelement)
```

### `Element` roll axis
Elements have a default roll angle with respect to its longitudinal axis of $\pi/2$, which corresponds to keeping the strong bending axis flat against the XY plane. If you wish to change this, you can change it by accessing the Ψ parameter:
```julia
element.Ψ = pi
```

### `Element` release
Elements can have partial DOF releases to decouple nodal displacements from element end displacements. This is the process of adding hinges to one or both ends of the beam. You can do this in the construction of an element through the optional argument `release`:
```julia
released_element = Element(pin_support, roller_support, ibeam_section; release = :freefixed)
```

By default, no releases are performed (i.e., `release = :fixedfixed`). You can choose between:
- `:freefixed` create a hinge in the beginning node of the element
- `:fixedfree` create a hinge in the ending node of the element
- `:freefree` create hinges on both ends of the element
- `:joist` create hinges on both ends of the element with the exception of torsional DOFs.

# Loads
Loads can be applied to nodes and elements.

## `NodeForce`
A `NodeForce` is defined on a node with a force vector:
```julia
load1 = NodeForce(free_nodes[1], [0., 0., -150.0])
```

## `NodeMoment`
A `NodeMoment` is defined on a node with a moment vector ($M_x, M_y, M_z$):
```julia
load2 = NodeForce(free_nodes[5], [40., 0., 0.])
```

## `LineLoad`
A `LineLoad` is defined on an element with a force vector, whose magnitude indicates the length-normalized force value, $\text{force}/\text{distance}$. E.g. a downwards load of $10\text{kN}/\text{m}$ is defined as (assuming we are working in kN, m):
```julia
snow_load = LineLoad(element, [0., 0., -10.])
```

## `PointLoad`
A `PointLoad` is defined on an element with a normalized position $0<x<1$ where the load is applied and the load value. E.g., a load of $20$ in the X axis direction applied at the quarter point of an element is defined as:
```julia
sideways_load = PointLoad(element, 0.25, [20.0, 0., 0.])
```

## `Model`
```julia
mutable struct Model{E,L} <: AbstractModel
    nodes::Vector{Node}
    elements::Vector{E}
    loads::Vector{L}
    nNodes::Int64
    nElements::Int64
    DOFs::Vector{Bool} #vector of DOFs
    nDOFs::Int64
    freeDOFs::Vector{Int64} #free DOF indices
    fixedDOFs::Vector{Int64}
    S::SparseMatrixCSC{Float64,Int64} # global stiffness
    P::Vector{Float64} # external loads
    Pf::Vector{Float64} # element end forces
    u::Vector{Float64} # nodal displacements
    reactions::Vector{Float64} # reaction forces
    compliance::Float64 #structural compliance
    tol::Float64
    processed::Bool
end
```

A model is assembled from a collection of nodes, elements, and loads:
```julia
nodes = [pin_support; roller_support; free_nodes]
elements = [element, released_element]
loads = [load1, load2, snow_load, sideways_load]

model = Model(nodes, elements, loads)
```

## Solving
The primary unknown field we are trying to find is `u`, the vector of all nodal DOFs (in order of assembly) in which equilibrium holds. We can find this via:
```julia
solve!(model)
```
(Note that in this nonsensical example, this will result in a singular error).

You can access the solved field via:
```julia
u = model.u
```

Or directly from the populated fields in the nodes:
```julia
node2_displacement = model.nodes[2].displacement
```

If a node has a restrained DOF, you can find its reaction from:
```julia
roller_reaction_forces = roller_support.reaction
```

You can also find the end forces acting on an element via:
```julia
element_forces = model.elements[2].forces
```
Which gives a vector: $F_{x1}, F_{y1}, F_{z1}, M_{x1}, M_{y1}, M_{z1}, F_{x2}, F_{y2}, F_{z2}, M_{x2}, M_{y2}, M_{z2}$ where $1$ is the starting node and $2$ is the ending node, with all values defined in the *local coordinate system* of the beam.

## New loads
If you have a new set of loads, directly get the corresponding displacement via:
```julia
u_new = solve!(model, new_loads)
```
Or replace the vector of loads associated with the model and solve in place via:
```julia
solve!(model, new_loads)
```

## Updating values
If you change a value, such as the position of a node, reprocess the fields before solving by:
```julia
#change 1
model.nodes[2].position .+= [5, 0, 0]
solve!(model;reprocess = true)
```

# Trusses
For truss structures, with only 3 translational DOFs per node, there are separate data structures for `TrussNode`, `TrussElement`, `TrussSection`, and `TrussModel`, which can be defined similarily as above except:

1. `TrussNode`s are constructed using only a length 3 vector of booleans if you are explicitly defining the DOF restrictions: `TrussNode(rand(3), [true, true, false])`.
2. `TrussSection`s only require the area and length, `TrussSection(A, E)`. You can use a regular `Section` to define a `TrussElement`.
3. `TrussElement`s do not have releases or roll angles. By definition they are equivalent to `Element(...; release = :freefree)`
4. Only `NodeForce`s can be applied as loads for `TrussModel`s.

## `TrussNode`
```julia
mutable struct TrussNode <: AbstractNode
    position::Vector{Float64}
    dof::Vector{Bool}
    nodeID::Int64
    globalID::Vector{Int64}
    reaction::Vector{Float64}
    displacement::Vector{Float64}
    id::Symbol
end
```

## `TrussElement`
```julia
mutable struct TrussElement <: AbstractElement
    section::Union{TrussSection,Section} #cross section
    nodeStart::TrussNode #start position
    nodeEnd::TrussNode #end position
    elementID::Int64
    globalID::Vector{Int64} #element global DOFs
    length::Float64 #length of element
    K::Matrix{Float64} # stiffness matrix in GCS
    R::Matrix{Float64} # transformation matrix
    forces::Vector{Float64} #elemental forces in LCS
    Ψ::Float64
    LCS::Vector{Vector{Float64}}
    id::Union{Symbol, Nothing} #optional identifier
end
```

## `TrussModel`
```julia
mutable struct TrussModel <: AbstractModel
    nodes::Vector{TrussNode}
    elements::Vector{TrussElement}
    loads::Vector{NodeForce}
    nNodes::Int64
    nElements::Int64
    DOFs::Vector{Bool} #vector of DOFs
    nDOFs::Int64
    freeDOFs::Vector{Int64} #free DOF indices
    fixedDOFs::Vector{Int64}
    S::SparseMatrixCSC{Float64,Int64} # global stiffness
    P::Vector{Float64} # external loads
    u::Vector{Float64} # nodal displacements
    reactions::Vector{Float64} # reaction forces
    compliance::Float64 #structural compliance
    tol::Float64
    processed::Bool
end
```