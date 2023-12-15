"""
    k_fixedfixed(element::Element)

[12 × 12] stiffness matrix for a beam element fully coupled at both ends.
"""
function k_fixedfixed(element::Element)
    E = element.section.E
    A = element.section.A
    L = element.length
    G = element.section.G
    Ix = element.section.Ix
    Iy = element.section.Iy
    J = element.section.J


    k = E / L^3 * [
        A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 12Ix 0 0 0 6L*Ix 0 -12Ix 0 0 0 6L*Ix;
        0 0 12Iy 0 -6L*Iy 0 0 0 -12Iy 0 -6L*Iy 0;
        0 0 0 G*J*L^2/E 0 0 0 0 0 -G*J*L^2/E 0 0;
        0 0 -6L*Iy 0 4L^2*Iy 0 0 0 6L*Iy 0 2L^2*Iy 0;
        0 6L*Ix 0 0 0 4L^2*Ix 0 -6L*Ix 0 0 0 2L^2*Ix;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -12Ix 0 0 0 -6L*Ix 0 12Ix 0 0 0 -6L*Ix;
        0 0 -12Iy 0 6L*Iy 0 0 0 12Iy 0 6L*Iy 0;
        0 0 0 -G*J*L^2/E 0 0 0 0 0 G*J*L^2/E 0 0;
        0 0 -6L*Iy 0 2L^2*Iy 0 0 0 6L*Iy 0 4L^2*Iy 0;
        0 6L*Ix 0 0 0 2L^2*Ix 0 -6L*Ix 0 0 0 4L^2*Ix
        ]

    return k
end

"""
    k_freefixed(element::Element)

[12 × 12] stiffness matrix for a beam element with rotational DOFs decoupled at the start node.
"""
function k_freefixed(element::Element)
    E = element.section.E
    A = element.section.A
    L = element.length
    Ix = element.section.Ix
    Iy = element.section.Iy

    k = E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 3Ix 0 0 0 0 0 -3Ix 0 0 0 3L*Ix;
        0 0 3Iy 0 0 0 0 0 -3Iy 0 -3L*Iy 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -3Ix 0 0 0 0 0 3Ix 0 0 0 -3L*Ix;
        0 0 -3Iy 0 0 0 0 0 3Iy 0 3L*Iy 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -3L*Iy 0 0 0 0 0 3L*Iy 0 3L^2*Iy 0;
        0 3L*Ix 0 0 0 0 0 -3L*Ix 0 0 0 3L^2*Ix    
    ]

    return k
end

"""
    k_fixedfree(element::Element)

[12 × 12] stiffness matrix for a beam element with rotational DOFs decoupled at the end node.
"""
function k_fixedfree(element::Element)
    E = element.section.E
    A = element.section.A
    L = element.length
    Ix = element.section.Ix
    Iy = element.section.Iy

    k = E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 3Ix 0 0 0 3L*Ix 0 -3Ix 0 0 0 0;
        0 0 3Iy 0 -3L*Iy 0 0 0 -3Iy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -3L*Iy 0 3L^2*Iy 0 0 0 3L*Iy 0 0 0;
        0 3L*Ix 0 0 0 3L^2*Ix 0 -3L*Ix 0 0 0 0 ;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -3Ix 0 0 0 -3L*Ix 0 3Ix 0 0 0 0;
        0 0 -3Iy 0 3L*Iy 0 0 0 3Iy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0    
    ]

    return k
end

"""
    k_freefree(element::Element)

[12 × 12] stiffness matrix for a beam element with full decoupled rotational DOFs.
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

"""
    k_joist(element::Element)

[12 × 12] stiffness matrix for a beam element with only torsional DOFs coupled to nodes.
"""
function k_joist(element::Element)
    E = element.section.E
    A = element.section.A
    L = element.length
    G = element.section.G
    J = element.section.J

    k = E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 G*J*L^2/E 0 0 0 0 0 -G*J*L^2/E 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 ;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 -G*J*L^2/E 0 0 0 0 0 G*J*L^2/E 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0    
    ]
end

const releases::Vector{Symbol} = [:fixedfixed, :freefixed, :fixedfree, :freefree, :joist, :truss]

const kDict = Dict(:fixedfixed => k_fixedfixed,
    :freefixed => k_freefixed,
    :fixedfree => k_fixedfree,
    :freefree => k_freefree,
    :truss => k_freefree,
    :joist => k_joist)

"""
    global_K!(element::Element)

Populate the element stiffness matrix `element.K` in GCS.
"""
function global_K!(element::Element)
    kfunction = kDict[element.release]
    element.K = Symmetric(element.R' * kfunction(element) * element.R)
end

function global_K!(element::TrussElement)
    k = element.section.E * element.section.A / element.length .* [1 -1; -1 1]
    element.K = Symmetric(element.R' * k * element.R)
end

"""
    global_K(element::Element)

Return the element stiffness matrix in GCS.
"""
function global_K(element::Element)
    kfunction = kDict[element.release]
    return Symmetric(element.R' * kfunction(element) * element.R)
end

function global_K(element::TrussElement)
    k = element.section.E * element.section.A / element.length .* [1 -1; -1 1]
    return Symmetric(element.R' * k * element.R)
end

"""
    local_K(element::Element)

Return the element stiffness matrix in LCS.
"""
function local_K(element::Element)
    kfunction = kDict[element.release]
    return kfunction(element)
end

function local_K(element::TrussElement)
    return element.section.E * element.section.A / element.length .* [1 -1; -1 1]
end
