abstract type AbstractNode end

"""
    Node(position::Vector{Float64}, dofs::Vector{Bool})

Instantiate a 6 DOF node with given position and fixities.

# Example
```julia-repl
julia> Node([4.3, 2.2, 10.4], [true, true, false, true, false, false])
Node([4.3, 2.2, 10.4], Bool[1, 1, 0, 1, 0, 0], #undef, #undef, #undef, nothing)
```
---------------------------------------------

    Node(position::Vector{Float64}, fixity::Symbol)

Instantiate a 6 DOF node with given position and common boundary type.
Available boundary conditions:
- :free
- :fixed
- :pinned
- :(x/y/z)free
- :(x/y/z)fixed
# Example
```julia-repl
julia> Node([4.3, 2.2, 10.4], :zfixed)
Node([4.3, 2.2, 10.4], Bool[1, 1, 0, 1, 1, 1], #undef, #undef, #undef, nothing)

```
"""
mutable struct Node <: AbstractNode
    position::Vector{Float64}
    dof::Vector{Bool}
    nodeID::Int64
    globalID::Vector{Int64}
    reaction::Vector{Float64}
    displacement::Vector{Float64}
    id::Union{Symbol, Nothing}

    function Node(position::Vector{Float64}, dofs::Vector{Bool})
        @assert length(position) == 3 && length(dofs) == 6 "Position vector must be in R³, DOFs must be length 6"

        node = new(position, dofs)

        node.id = nothing
        return node
    end

    function Node(position::Vector{Float64}, fixity::Symbol)

        @assert length(position) == 3 "Position vector must be in R³"

        dofs = copy(fixDict[fixity])

        node = new(position, dofs)

        node.id = nothing
        return node
    end
end

"""
    TrussNode(position::Vector{Float64}, dofs::Vector{Bool})

Instantiate a 3 DOF node with given position and fixities.

# Example
```julia-repl
julia> TrussNode([1., 1., 56.], [false, true, true])
TrussNode([1.0, 1.0, 56.0], Bool[0, 1, 1], #undef, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], nothing)
```
---------------------------------------------

    TrussNode(position::Vector{Float64}, fixity::Symbol)

Instantiate a 3 DOF node with given position and common boundary type.
Available boundary conditions:
- :free
- :pinned
- :(x/y/z)free
- :(x/y/z)fixed
# Example
```julia-repl
julia> TrussNode([1., 1., 56.], :pinned)
TrussNode([1.0, 1.0, 56.0], Bool[0, 0, 0], #undef, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], nothing)

```
"""
mutable struct TrussNode <: AbstractNode
    position::Vector{Float64}
    dof::Vector{Bool}
    nodeID::Int64
    globalID::Vector{Int64}
    reaction::Vector{Float64}
    displacement::Vector{Float64}
    id::Union{Symbol, Nothing}

    function TrussNode(position::Vector{Float64}, dofs::Vector{Bool})
        
        @assert length(position) == length(dofs) == 3  "Position and dof vector must be in R³"

        node = new(position, dofs)

        node.displacement = zeros(6)
        node.reaction = zeros(6)

        node.id = nothing
        return node
    end

    function TrussNode(position::Vector{Float64}, fixity::Symbol)
        
        @assert length(position) == 3 "Position vector must be in R³"

        dofs = copy(fixDict[fixity][1:3])

        node = new(position, dofs)

        node.displacement = zeros(3)
        node.reaction = zeros(3)

        node.id = nothing
        return node
    end
end

"""
Common fixity types
"""
const fixDict = Dict(:fixed => [false, false, false, false, false, false],
    :free => [true, true, true, true, true, true],
    :xfixed => [false, true, true, true, true, true],
    :yfixed => [true, false, true, true, true, true],
    :zfixed => [true, true, false, true, true, true],
    :xfree => [true, false, false, false, false, false],
    :yfree => [false, true, false, false, false, false],
    :zfree => [false, false, true, false, false, false],
    :pinned => [false, false, false, true, true, true])

"""
Inactive DOF w/r/t plane
"""
const planeDict = Dict(:XY => [3, 4, 5],
    :YZ => [1, 5, 6],
    :ZX => [2, 4, 6])

