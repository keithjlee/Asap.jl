# Code related to stiffness matrices

"""
Truss LCS K
"""
function K_truss(element::Element)
    return element.E * element.A / element.length .* [1 -1; -1 1]
end

"""
Frame LCS 2d K
"""
function K2d_frame(element::Element)
    # stiffness matrix
    k = zeros(6,6)

    p = [element.A / element.length, 0, 0, -element.A/element.length, 0, 0]
    v = [0, 12 * element.Iz/element.length^3, 6 * element.Iz/element.length^2, 0, -12element.Iz/element.length^3, 6element.Iz/element.length^2]
    m = [0, 6element.Iz/element.length^2, 4element.Iz/element.length, 0, -6element.Iz/element.length^2, 2element.Iz/element.length]

    k[1,:] = p
    k[2,:] = v
    k[3,:] = m
    k[4,:] = -p
    k[5,:] = -v
    k[6,:] = [0, 6element.Iz/element.length^2, 2element.Iz/element.length, 0, -6element.Iz/element.length^2, 4element.Iz/element.length]

    k .*= element.E

    return k
end

"""
3d frame K (dis is nasty)
"""
function K3d_frame(element::Element)
    return element.E / element.length^3 * [
        element.A*element.length^2 0 0 0 0 0 -element.A*element.length^2 0 0 0 0 0;
        0 12element.Iz 0 0 0 6element.length*element.Iz 0 -12element.Iz 0 0 0 6element.length*element.Iz;
        0 0 12element.Iy 0 -6element.length*element.Iy 0 0 0 -12element.Iy 0 -6element.length*element.Iy 0;
        0 0 0 element.G*element.J*element.length^2/element.E 0 0 0 0 0 -element.G*element.J*element.length^2/element.E 0 0;
        0 0 -6element.length*element.Iy 0 4element.length^2*element.Iy 0 0 0 6element.length*element.Iy 0 2element.length^2*element.Iy 0;
        0 6element.length*element.Iz 0 0 0 4element.length^2*element.Iz 0 -6element.length*element.Iz 0 0 0 2element.length^2*element.Iz;
        -element.A*element.length^2 0 0 0 0 0 element.A*element.length^2 0 0 0 0 0;
        0 -12element.Iz 0 0 0 -6element.length*element.Iz 0 12element.Iz 0 0 0 -6element.length*element.Iz;
        0 0 -12element.Iy 0 6element.length*element.Iy 0 0 0 12element.Iy 0 6element.length*element.Iy 0;
        0 0 0 -element.G*element.J*element.length^2/element.E 0 0 0 0 0 element.G*element.J*element.length^2/element.E 0 0;
        0 0 -6element.length*element.Iy 0 2element.length^2*element.Iy 0 0 0 6element.length*element.Iy 0 4element.length^2*element.Iy 0;
        0 6element.length*element.Iz 0 0 0 2element.length^2*element.Iz 0 -6element.length*element.Iz 0 0 0 4element.length^2*element.Iz
        ]
end

"""
matches condition to appropriate stiffness matrix
"""
Kdict = Dict(
    (2, :truss) => K_truss,
    (3, :truss) => K_truss,
    (2, :frame) => K2d_frame,
    (3, :frame) => K3d_frame
)

"""
Creates elemental stiffness matrix in GCS
"""
function k_elemental!(element::Element, dims::Int; return_k = false, tol = 1e-3)

    rFunction = Rdict[(dims, element.type)] # rotation matrix generator
    kFunction = Kdict[(dims, element.type)] # stiffness matrix generator

    rFunction(element) # add R to element.R
    K = kFunction(element) # LCS stiffness matrix
    element.k = Symmetric(element.R' * K * element.R) #GCS stiffness matrix

    if return_k
        return element.k
    end
end

"""
Creates global stiffness matrix K
"""
function K(elements::Array{Element}, nDOFS::Int)
    K = spzeros(nDOFS, nDOFS)

    for element in elements
        K[element.dofIndex, element.dofIndex] .+= element.k
    end

    return K
end