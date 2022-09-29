# Custom indexing
function Base.getindex(elements::Vector{Element}, i::Symbol)
    return [element for element in elements if element.id == i]
end

function Base.findall(elements::Vector{Element}, i::Symbol)
    return findall([element.id == i for element in elements])
end

function Base.getindex(nodes::Vector{Node}, i::Symbol)
    return [node for node in nodes if node.id == i]
end

function Base.findall(nodes::Vector{Node}, i::Symbol)
    return findall([node.id == i for node in nodes])
end

function Base.getindex(loads::Vector{Load}, i::Symbol)
    return [load for load in loads if load.id == i]
end

function Base.findall(loads::Vector{Load}, i::Symbol)
    return findall([load.id == i for load in loads])
end

function Base.getindex(structure::Structure, i::Symbol)
    return structure.nodes[i], structure.elements[i], structure.loads[i]
end