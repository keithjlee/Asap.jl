function Base.getindex(loads::Vector{Load}, i::Symbol)
    return [load for load in loads if load.id == i]
end

function Base.findall(loads::Vector{Load}, i::Symbol)
    return findall([load.id == i for load in loads])
end