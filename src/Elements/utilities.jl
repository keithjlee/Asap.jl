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
Local coordinate system of element with roll angle
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
displacement function for the transverse translation and in-plane rotation for a GIVEN PLANE IN THE LCS OF AN ELEMENT:

u_xy= [ v₁
        θ₁
        v₂
        θ₂ ]

IE for the local XY plane:

v₁, v₂ are the start and end displacements in the local Y direction
θ₁, θ₂ are the start and end rotations in the **local Z** direction (ie rotation in plane of local XY)

Gives:

v_y(x) = N × u_xy (translational displacement in local Y at point x)

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

"""
Second derivative of (assumed cubic) shape function N, for use in calculating internal moments:

M(x) = -EI d²v/dx² = -EI ⋅ (d²N/dx² × u_xy)

Note that this is essentially useless as it reduces a cubic function into a linear one and thus cannot accurately decribe the moment across an element unless highly discretized.

As such the third derivative (d³N/dx³) to determine internal Shear is also effectively useless as it will give a constant value across the entire section
"""
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
Axial displacement function: linear interpolation between start and end displacements
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

"""
    displacements(element::Element; n::Integer = 20)

Get the [3 × n] matrix where each column represents the local [x,y,z] displacement of the element from end forces
"""
function displacements(element::Element; n::Integer = 20)

    # base properties
    ulocal = element.R * [element.nodeStart.displacement; element.nodeEnd.displacement]
    # d = [element.nodeStart.dof; element.nodeEnd.dof]
    # ulocal[.!d] .= 0
    L = element.length

    # extracting relevant nodal DOFs
    uX = ulocal[[1, 7]]
    uY = ulocal[[2, 6, 8, 12]]
    uZ = ulocal[[3, 5, 9, 11]]

    # discretizing length of element
    xrange = range(0, L, n)

    [@inbounds hcat([Naxial(x, L) * uX for x in xrange]...);
        @inbounds hcat([N(x, L) * uY for x in xrange]...);
        @inbounds hcat([N(x, L) * uZ for x in xrange]...)]
end


"""
    displacedshape(element::Element; factor::Real = 1, n::Integer = 20)

Generate a [3 × n] matrix whose columns describe the displaced shape of an element with a displacement scale factor of `factor`
"""
function displacedshape(element::Element; factor::Real = 1, n::Integer = 20)
    Δ = displacements(element; n = n)

    xlocal = first(element.LCS)

    init = element.nodeStart.position

    inc = range(0, element.length, n)

    hcat([init .+ xlocal .* i .+ factor .* sum(disp .* element.LCS) for (i, disp) in zip(inc, eachcol(Δ))]...)
end