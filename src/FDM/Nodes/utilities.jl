# function Base.-(nodeEnd::FDMnode, nodeStart::FDMnode)
#     nodeEnd.position - nodeStart.position
# end

"""
Custom indexing based on IDs of structs
"""
function Base.getindex(nodes::Vector{FDMnode}, i::Symbol)
    return [node for node in nodes if node.id == i]
end

function Base.findall(nodes::Vector{FDMnode}, i::Symbol)
    return findall([x.id == i for x in nodes])
end