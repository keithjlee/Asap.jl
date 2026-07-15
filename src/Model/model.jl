"""
    Model{T<:Real}

A structural model: the complete *definition* of a structure — nodes,
elements (frame and truss mixed freely), loads, and elastic supports.

The model is the first of three separated layers:

1. **`Model`** — user definition data (this type). Mutable containers so
   models can be built and edited incrementally.
2. **`AnalysisCache`** — analysis structure derived from the definition by
   [`process!`](@ref): DOF partition, sparsity pattern, buffers, cached
   factorization. Rebuilt when topology changes.
3. **`LinearResults`** — outputs of [`solve!`](@ref): displacements,
   reactions, element forces, queried through accessor functions.

# Fields
- `nodes::Vector{Node{T}}`
- `elements::Vector{AbstractElement{T}}`: heterogeneous by design — truss
  and frame elements coexist; assembly groups them by concrete type behind
  function barriers so hot loops stay type-stable
- `loads::Vector{AbstractLoad{T}}`
- `springs::Vector{NodalSpring{T}}`: elastic supports (applicative — nodes
  don't know about their springs)
- `cache`: the `AnalysisCache` after processing (`nothing` before)
- `results`: the most recent `LinearResults` (`nothing` before solving)

# Constructor
    Model(nodes, elements, loads = AbstractLoad{T}[]; springs = NodalSpring{T}[])

# Typical use
```julia-repl
julia> model = Model(nodes, elements, loads)

julia> solve!(model)

julia> displacement(model.results, nodes[end])
```
"""
mutable struct Model{T<:Real}
    nodes::Vector{Node{T}}
    elements::Vector{AbstractElement{T}}
    loads::Vector{AbstractLoad{T}}
    springs::Vector{NodalSpring{T}}
    cache::Any        # ::Union{Nothing, AnalysisCache{T}} — Any avoids a mutual-recursion forward reference
    results::Any      # ::Union{Nothing, LinearResults{T}}

    function Model(nodes::Vector{Node{T}}, elements::Vector{<:AbstractElement{T}},
        loads::Vector{<:AbstractLoad{T}}=AbstractLoad{T}[];
        springs::Vector{NodalSpring{T}}=NodalSpring{T}[]) where {T}
        @assert !isempty(nodes) "a model needs nodes"
        @assert !isempty(elements) "a model needs elements"
        return new{T}(nodes, collect(AbstractElement{T}, elements),
            collect(AbstractLoad{T}, loads), springs, nothing, nothing)
    end
end

"""
    node_positions(model) -> Matrix{T}

`n × 3` matrix of nodal coordinates (one row per node).
"""
node_positions(model::Model{T}) where {T} =
    permutedims(reduce(hcat, (collect(n.position) for n in model.nodes)))

"""
    connectivity(model) -> SparseMatrixCSC

`n_elements × n_nodes` signed incidence matrix: −1 at each element's start
node, +1 at its end node — so `C * X` gives element vectors for nodal
coordinate columns `X`.
"""
function connectivity(model::Model)
    I_ = Int[]
    J_ = Int[]
    V_ = Int[]
    for (i, el) in enumerate(model.elements)
        append!(I_, (i, i))
        append!(J_, (el.nodeStart.index, el.nodeEnd.index))
        append!(V_, (-1, 1))
    end
    return sparse(I_, J_, V_, length(model.elements), length(model.nodes))
end

"""
    volume(model) -> T

Total material volume: Σ A·L over elements whose sections expose a
geometric area (`Section`); `RigiditySection`s have no geometric area and
contribute zero (they represent effective rigidities, not geometry).
"""
volume(model::Model{T}) where {T} =
    sum(el -> _section_area(el.section) * Base.length(el), model.elements; init=zero(T))

_section_area(s::Section) = s.A
_section_area(::AbstractSection{T}) where {T} = zero(T)

# indexing sugar carried over from the legacy API: query nodes/elements by id
Base.getindex(nodes::Vector{<:Node}, id::Symbol) = [n for n in nodes if n.id == id]
Base.findall(nodes::Vector{<:Node}, id::Symbol) = [i for (i, n) in enumerate(nodes) if n.id == id]
Base.getindex(els::Vector{<:AbstractElement}, id::Symbol) = [e for e in els if e.id == id]
Base.findall(els::Vector{<:AbstractElement}, id::Symbol) = [i for (i, e) in enumerate(els) if e.id == id]

"""
    planarize!(model, plane = :XY)

Constrain the model to 2D behavior: planarize all nodes (fix out-of-plane
DOFs) and, for the `:XY` plane, zero every frame element's roll angle `rollangle`
so local bending axes align with the working plane (legacy behavior,
preserved).
"""
function planarize!(model::Model, plane::Symbol=:XY)
    planarize!(model.nodes, plane)
    if plane == :XY
        for el in model.elements
            el isa FrameElement && (el.rollangle = zero(el.rollangle))
        end
    end
    return model
end
