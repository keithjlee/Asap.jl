"""
    k_fixedfixed(element::Element)

[12 × 12] stiffness matrix for a beam element fully coupled at both ends.
"""
function k_fixedfixed(E, A, L, G, Ix, Iy, J)

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
function k_freefixed(E, A, L, Ix, Iy)

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
function k_fixedfree(E, A, L, Ix, Iy)
    
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
function k_freefree(E, A, L)

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
function k_joist(E, A, L, G, J)

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

    return k
end

"""
    local_K(element::Element)

Return the element stiffness matrix in LCS.
"""
local_K(element::Element{FixedFixed}) = k_fixedfixed(element.section.E, element.section.A, element.length, element.section.G, element.section.Ix, element.section.Iy, element.section.J)
local_K(element::Element{FixedFree}) = k_fixedfree(element.section.E, element.section.A, element.length, element.section.Ix, element.section.Iy)
local_K(element::Element{FreeFixed}) = k_freefixed(element.section.E, element.section.A, element.length, element.section.Ix, element.section.Iy)
local_K(element::Element{FreeFree}) = k_freefree(element.section.E, element.section.A, element.length)
local_K(element::Element{Joist}) = k_joist(element.section.E, element.section.A, element.length, element.section.G, element.section.J)

local_K(element::TrussElement) = element.section.E * element.section.A / element.length .* [1 -1; -1 1]

"""
    global_K!(element::Element)

Populate the element stiffness matrix `element.K` in GCS.
"""
function global_K!(element::Element)
    element.K = element.R' * local_K(element) * element.R
end

function global_K!(element::TrussElement)
    element.K = element.R' * local_K(element) * element.R
end

"""
    global_K(element::Element)

Return the element stiffness matrix in GCS.
"""
function global_K(element::Element)
    return element.R' * local_K(element) * element.R
end

function global_K(element::TrussElement)
    return element.R' * local_K(element) * element.R
end