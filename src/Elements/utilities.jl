"""
Extract all elements with a given ID
elements = Vector{Element}()
elements[:outerEdges]
"""
function Base.getindex(elements::Vector{<:AbstractElement}, i::Symbol)
    return [element for element in elements if element.id == i]
end

function Base.findall(elements::Vector{<:AbstractElement}, i::Symbol)
    return findall([element.id == i for element in elements])
end

function Base.getindex(elements::Vector{TrussElement}, i::Symbol)
    return [element for element in elements if element.id == i]
end

function Base.findall(elements::Vector{TrussElement}, i::Symbol)
    return findall([element.id == i for element in elements])
end

Base.length(element::AbstractElement) = norm(element.nodeEnd.position .- element.nodeStart.position)

function length!(element::AbstractElement)
    element.length = length(element)
end

"""
    release!(element::Element, release::Symbol)

Change the release condition of an element.
Available releases:
- :fixedfixed (default)
- :fixedfree
- :freefixed
- :freefree
"""
function release!(element::Element, release::Symbol)
    @assert in(release, releases) "Release not recognized; choose from: :fixedfixed, :freefixed, :fixedfree, :freefree"

    element.release = release
    element.Q = zeros(12)
    global_K!(element)
end


"""
    local_x(element::AbstractElement; unit = true)

Get the local x vector of an element: element.nodeEnd.position - element.nodeStart.position.
`unit = true` gives the normalized vector.
"""
function local_x(element::AbstractElement; unit = true)
    x = element.nodeEnd.position .- element.nodeStart.position
    unit ? normalize(x) : x
end

"""
    lcs(element::AbstractElement, Ψ::Float64; tol = 1e-6)

Get the local coordinate system unit vectors of a given element and pitch angle Ψ: [local_x, local_y, local_z]
"""
function lcs(element::AbstractElement, Ψ::Float64; tol = 1e-6)

    # local x vector
    xvec = local_x(element)
    
    if norm(cross(xvec, globalY)) < tol
        CYx = xvec[2] #cosine to global Y axis
        xvec = CYx * globalY
        yvec = -CYx * globalX * cos(Ψ) + sin(Ψ) * globalZ
        zvec = CYx * globalX * sin(Ψ) + cos(Ψ) * globalZ 
    else
        zbar = normalize(cross(xvec, [0, 1, 0]))
        ybar = normalize(cross(zbar, xvec))

        yvec = cos(Ψ) * ybar + sin(Ψ) * zbar
        zvec = -sin(Ψ) * ybar + cos(Ψ) * zbar
    end

    return [xvec, yvec, zvec]
end

"""
    lcs(element::AbstractElement, Ψ::Float64; tol = 1e-6)

Populate local coordinate system unit vectors of a given element and pitch angle Ψ: [local_x, local_y, local_z]
"""
function lcs!(element::AbstractElement, Ψ::Float64; tol = 0.001)

    # local x vector
    xvec = local_x(element)
    
    if norm(cross(xvec, globalY)) < tol
        CYx = xvec[2] #cosine to global Y axis
        xvec = CYx * globalY
        yvec = -CYx * globalX * cos(Ψ) + sin(Ψ) * globalZ
        zvec = CYx * globalX * sin(Ψ) + cos(Ψ) * globalZ 
    else
        zbar = normalize(cross(xvec, [0, 1, 0]))
        ybar = normalize(cross(zbar, xvec))

        yvec = cos(Ψ) * ybar + sin(Ψ) * zbar
        zvec = -sin(Ψ) * ybar + cos(Ψ) * zbar
    end

    element.LCS = [xvec, yvec, zvec]
end


"""
    lcs(xin::Vector{Float64}, Ψ::Float64; tol = 1e-4)

Get the local coordinate system unit vectors of a given local x axis and pitch angle Ψ: [local_x, local_y, local_z]
"""
function lcs(xin::Vector{Float64}, Ψ::Float64; tol = 1e-4)

    xvec = normalize(copy(xin))
    # local x vector
    if norm(cross(xvec, globalY)) < tol
        CYx = xvec[2] #cosine to global Y axis
        xvec = CYx * globalY
        yvec = -CYx * globalX * cos(Ψ) + sin(Ψ) * globalZ
        zvec = CYx * globalX * sin(Ψ) + cos(Ψ) * globalZ 
    else
        zbar = normalize(cross(xvec, [0, 1, 0]))
        ybar = normalize(cross(zbar, xvec))

        yvec = cos(Ψ) * ybar + sin(Ψ) * zbar
        zvec = -sin(Ψ) * ybar + cos(Ψ) * zbar
    end

    return [xvec, yvec, zvec]
end

"""
    endpoints(element::AbstractElement)

Extract the start and end points as two vectors.
"""
function endpoints(element::AbstractElement)
    return [element.nodeStart.position, element.nodeEnd.position]
end

"""
    midpoint(element::AbstractElement)

Extract the centerpoint of an element as a vector.
"""
function midpoint(element::AbstractElement)
    return (element.nodeStart.position .+ element.nodeEnd.position) ./ 2
end

"""
    axial_force(element::TrussElement)

Extract the axial force of an element
"""
function axial_force(element::TrussElement)
    element.forces[2]
end

"""
    axial_force(element::Element)
    
Extract the axial force of an element
"""
function axial_force(element::Element)
    element.forces[7]
end