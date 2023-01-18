# global axes
const globalX = [1., 0., 0.]
const globalY = [0., 1., 0.]
const globalZ = [0., 0., 1.]

"""
k for frame element
"""
function k_fixedfixed(element::Element)
    E = element.section.E
    A = element.section.A
    L = element.length
    G = element.section.G
    Izz = element.section.Izz
    Iyy = element.section.Iyy
    J = element.section.J


    k = E / L^3 * [
        A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 12Izz 0 0 0 6L*Izz 0 -12Izz 0 0 0 6L*Izz;
        0 0 12Iyy 0 -6L*Iyy 0 0 0 -12Iyy 0 -6L*Iyy 0;
        0 0 0 G*J*L^2/E 0 0 0 0 0 -G*J*L^2/E 0 0;
        0 0 -6L*Iyy 0 4L^2*Iyy 0 0 0 6L*Iyy 0 2L^2*Iyy 0;
        0 6L*Izz 0 0 0 4L^2*Izz 0 -6L*Izz 0 0 0 2L^2*Izz;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -12Izz 0 0 0 -6L*Izz 0 12Izz 0 0 0 -6L*Izz;
        0 0 -12Iyy 0 6L*Iyy 0 0 0 12Iyy 0 6L*Iyy 0;
        0 0 0 -G*J*L^2/E 0 0 0 0 0 G*J*L^2/E 0 0;
        0 0 -6L*Iyy 0 2L^2*Iyy 0 0 0 6L*Iyy 0 4L^2*Iyy 0;
        0 6L*Izz 0 0 0 2L^2*Izz 0 -6L*Izz 0 0 0 4L^2*Izz
        ]

    return k
end

"""
Hinge-Fixed
"""
function k_freefixed(element::Element)
    E = element.section.E
    A = element.section.A
    L = element.length
    Izz = element.section.Izz
    Iyy = element.section.Iyy

    k = E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 3Izz 0 0 0 0 0 -3Izz 0 0 0 3L*Izz;
        0 0 3Iyy 0 0 0 0 0 -3Iyy 0 -3L*Iyy 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -3Izz 0 0 0 0 0 3Izz 0 0 0 -3L*Izz;
        0 0 -3Iyy 0 0 0 0 0 3Iyy 0 3L*Iyy 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -3L*Iyy 0 0 0 0 0 3L*Iyy 0 3L^2*Iyy 0;
        0 3L*Izz 0 0 0 0 0 -3L*Izz 0 0 0 3L^2*Izz    
    ]

    return k
end

"""
Fixed-Hinge
"""
function k_fixedfree(element::Element)
    E = element.section.E
    A = element.section.A
    L = element.length
    Izz = element.section.Izz
    Iyy = element.section.Iyy

    k = E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 3Izz 0 0 0 3L*Izz 0 -3Izz 0 0 0 0;
        0 0 3Iyy 0 -3L*Iyy 0 0 0 -3Iyy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -3L*Iyy 0 3L^2*Iyy 0 0 0 3L*Iyy 0 0 0;
        0 3L*Izz 0 0 0 3L^2*Izz 0 -3L*Izz 0 0 0 0 ;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -3Izz 0 0 0 -3L*Izz 0 3Izz 0 0 0 0;
        0 0 -3Iyy 0 3L*Iyy 0 0 0 3Iyy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0    
    ]

    return k
end

"""
Hinge-Hinge
"""
function k_freefree(element::Element)
    E = element.section.E
    A = element.section.A
    L = element.length

    k = zeros(12, 12)
    k[1,1] = 1
    k[1,7] = -1
    k[7,1] = -1
    k[7,7] = 1

    return E * A / L .* k
end

releases = [:fixedfixed, :freefixed, :fixedfree, :freefree]

const kDict = Dict(:fixedfixed => k_fixedfixed,
    :freefixed => k_freefixed,
    :fixedfree => k_fixedfree,
    :freefree => k_freefree,
    :truss => k_freefree)

"""
make global k
"""
function makeK!(element::Element)
    kfunction = kDict[element.release]
    element.K = Symmetric(element.R' * kfunction(element) * element.R)
end

"""
make global k
"""
function makeK!(element::TrussElement)
    k = element.section.E * element.section.A / element.length .* [1 -1; -1 1]
    element.K = Symmetric(element.R' * k * element.R)
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
Transformation matrix
"""
function R(element::Union{Element, GeometricElement}; tol = 1e-4)
    xvec = normalize(element.posEnd .- element.posStart) # local x vector
    CXx, CYx, CZx = xvec # local x cosines

    cΨ = cos(element.Ψ)
    sΨ = sin(element.Ψ)


    if norm(cross(xvec, globalY)) < tol #special case for horizontal members aligned with global Y
        Λ = [0. CYx 0.;
            -CYx*cΨ 0 sΨ;
            CYx*sΨ 0 cΨ]
    else # all other
        b1 = (-CXx * CYx * cΨ - CZx * sΨ) / sqrt(CXx^2 + CZx^2)
        b2 = sqrt(CXx^2 + CZx^2) * cΨ
        b3 = (-CYx * CZx * cΨ + CXx * sΨ) / sqrt(CXx^2 + CZx^2)

        c1 = (CXx * CYx * sΨ - CZx * cΨ) / sqrt(CXx^2 + CZx^2)
        c2 = -sqrt(CXx^2 + CZx^2) * sΨ
        c3 = (CYx * CZx * sΨ + CXx * cΨ) / sqrt(CXx^2 + CZx^2)

        Λ = [CXx CYx CZx; 
            b1 b2 b3; 
            c1 c2 c3]
    end
    
    R = [Λ zeros(3,9); zeros(3,3) Λ zeros(3,6); zeros(3,6) Λ zeros(3,3); zeros(3,9) Λ]

    return R
end

"""
Transformation matrix
"""
function R(element::TrussElement)
    Cx, Cy, Cz = (element.posEnd .- element.posStart) ./ element.length
    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]
end


"""
Local coordinate system of element
"""
function lcs(element::Element, Ψ; tol = 0.001)

    # local x vector
    xvec = normalize(element.posEnd .- element.posStart)
    
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
displacement function
"""
function N(x::Float64, L::Float64)
    n1 = 1 - 3(x/L)^2 + 2(x/L)^3
    n2 = x*(1 - x/L)^2
    n3 = 3(x/L)^2 - 2(x/L)^3
    n4 = x^2/L * (-1 + x/L)

    return [n1 n2 n3 n4]
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