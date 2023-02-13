"""
K: Fixed-Fixed
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
K: Hinge-Fixed
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
K: Fixed-Hinge
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
K: Hinge-Hinge
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

const releases::Vector{Symbol} = [:fixedfixed, :freefixed, :fixedfree, :freefree]

const kDict = Dict(:fixedfixed => k_fixedfixed,
    :freefixed => k_freefixed,
    :fixedfree => k_fixedfree,
    :freefree => k_freefree,
    :truss => k_freefree)

"""
make elemental stiffness matrix in GCS
"""
function makeK!(element::Element)
    kfunction = kDict[element.release]
    element.K = Symmetric(element.R' * kfunction(element) * element.R)
end

function makeK(element::Element)
    kfunction = kDict[element.release]
    return Symmetric(element.R' * kfunction(element) * element.R)
end

function localK(element::Element)
    kfunction = kDict[element.release]
    return kfunction(element)
end

"""
make elemental stiffness matrix in GCS
"""
function makeK!(element::TrussElement)
    k = element.section.E * element.section.A / element.length .* [1 -1; -1 1]
    element.K = Symmetric(element.R' * k * element.R)
end

function localK(element::TrussElement)
    return element.section.E * element.section.A / element.length .* [1 -1; -1 1]
end
