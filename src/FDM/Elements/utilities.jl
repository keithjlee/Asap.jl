"""
    local_x(element; unit = true) -> SVector{3}

The element's direction vector at the CURRENT node positions (unit by
default).
"""
function local_x(element::FDMelement; unit=true)
    vec = element.pEnd.position - element.pStart.position
    unit ? normalize(vec) : vec
end

"member length at the current node positions [length]"
Base.length(element::FDMelement) = norm(local_x(element; unit=false))

"""
    force(element::FDMelement) -> T

Member axial force at the current geometry: `q · L` [force] (tension +
for positive force densities).
"""
force(element::FDMelement) = Base.length(element) * element.q

function Base.getindex(elements::Vector{<:FDMelement}, i::Symbol)
    return [element for element in elements if element.id == i]
end

function Base.findall(elements::Vector{<:FDMelement}, i::Symbol)
    return findall([x.id == i for x in elements])
end
