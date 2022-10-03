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
3d frame K
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

# creates the elemental stiffness matrix in global coordinate system
function k_elemental!(element::Element, dims::Int; return_k = false, tol = 1e-3)
    type = element.type
    if type == :truss
        if dims == 3
            R3d_truss!(element) # rotation matrix
            K = K_truss(element) # LCS stiffness matrix
            element.k = Symmetric(element.R' * K * element.R)
        elseif dims == 2
            R2d_truss!(element) # rotation matrix
            K = K_truss(element) # LCS stiffness matrix
            element.k = Symmetric(element.R' * K * element.R)
        else
            error("dims must be 2 or 3")
        end
    elseif type == :frame
        if dims == 3
            R3d_frame!(element) # rotation matrix
            K = K3d_frame(element) # LCS stiffness matrix
            element.k = Symmetric(element.R' * K * element.R)
        elseif dims == 2
            R2d_frame!(element) # rotation matrix
            K = K2d_frame(element) # LCS stiffness matrix
            element.k = Symmetric(element.R' * K * element.R)
        else
            error("dims must be 2 or 3")
        end
    else
        error("type must be :truss or :frame")
    end

    if return_k
        return element.k
    end
end

function K(elements::Array{Element}, nDOFS::Int)
    K = spzeros(nDOFS, nDOFS)

    for element in elements
        K[element.dofIndex, element.dofIndex] .+= element.k
    end

    return K
end