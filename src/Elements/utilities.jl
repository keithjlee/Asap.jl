"""
Extract all elements with a given ID
elements = Vector{Element}()
elements[:outerEdges]
"""
function Base.getindex(elements::Vector{Element}, i::Symbol)
    return [element for element in elements if element.id == i]
end

function Base.findall(elements::Vector{Element}, i::Symbol)
    return findall([element.id == i for element in elements])
end

function Base.getindex(elements::Vector{TrussElement}, i::Symbol)
    return [element for element in elements if element.id == i]
end

function Base.findall(elements::Vector{TrussElement}, i::Symbol)
    return findall([element.id == i for element in elements])
end


"""
change the release of an element
"""
function release!(element::Element, release::Symbol)
    if !in(release, releases)
        error("Release not recognized; choose from: :fixedfixed, :freefixed, :fixedfree, :freefree")
    end

    element.release = release
    element.Q = zeros(12)
    makeK!(element)
end


"""
Get the local x vector of an element
Defaults to normalized unit vector
"""
function localx(element::AbstractElement; unit = true)
    x = element.nodeEnd.position .- element.nodeStart.position
    unit ? normalize(x) : x
end

"""
Local coordinate system of element
"""
function lcs(element::AbstractElement, Ψ::Float64; tol = 0.001)

    # local x vector
    xvec = localx(element)
    
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
Local coordinate system of a directional vector and roll angle
"""
function lcs(xin::Vector{Float64}, Ψ::Float64; tol = 0.001)

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
Extract start and end positions of element
"""
function endPositions(element::AbstractElement)
    return [element.nodeStart.position, element.nodeEnd.position]
end


"""
displacement function
"""
function N(x::Float64, L::Float64)
    n1 = 1 - 3(x/L)^2 + 2(x/L)^3
    n2 = x*(1 - x/L)^2
    n3 = 3(x/L)^2 - 2(x/L)^3
    n4 = x^2/L * (-1 + x/L)

    # n1 = 1 + 1 / L^3 * (-3L * x + 2x^3)
    # n2 = x + 1 / L^2 * (-2L * x^2 + x^3)
    # n3 = 1 / L^3 * (3L * x^2 - 2x^3)
    # n4 = 1 / L^2 * (-L * x^2 + x^3)

    return [n1 n2 n3 n4]
end

function dN2(x::Float64, L::Float64)
    # n1 = 1 - 3(x/L)^2 + 2(x/L)^3
    # n2 = x*(1 - x/L)^2
    # n3 = 3(x/L)^2 - 2(x/L)^3
    # n4 = x^2/L * (-1 + x/L)

    n1 = 1 / L^3 * (12x - 6L)
    n2 = 1 / L^2 * (6x - 4L)
    n3 = 1 / L^3 * (6L - 12x)
    n4 = 1 / L^2 * (-2L + 6x)

    return [n1 n2 n3 n4]
end

"""
Axial displacement function
"""
function Naxial(x::Float64, L::Float64)
    n1 = 1 - x/L
    n2 = x / L

    return [n1 n2]
end

"""
stress function
"""
function B(y::Float64, x::Float64, L::Float64)
    b1 = 6 * (-1 + 2 * x / L)
    b2 = 2L * (-2 + 3 * x / L)
    b3 = 6 * (1 - 2 * x / L)
    b4 = 2L * (-1 + 3 * x / L)

    return -y/L^2 .* [b1 b2 b3 b4]
end