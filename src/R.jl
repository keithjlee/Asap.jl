# Rotation matrices

"""
2d truss rotation matrix
"""
function R2d_truss(element::Element)
    Cx, Cy = (element.posEnd .- element.posStart) ./ element.length
    return [Cx Cy 0 0 ; 0 0 Cx Cy]
end

"""
2d truss rotation matrix
"""
function R2d_truss!(element::Element)
    element.R = R2d_truss(element)
end

"""
3d truss rotation matrix
"""
function R3d_truss(element::Element)
    Cx, Cy, Cz = (element.posEnd .- element.posStart) ./ element.length
    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]
end

function R3d_truss!(element::Element)
    element.R = R3d_truss(element)
end

"""
2d frame rotation matrix
"""
function R2d_frame(element::Element)
    Cx, Cy = (element.posEnd .- element.posStart) ./ element.length

    # rotation matrix
    R = [Cx Cy 0 0 0 0;
        -Cy Cx 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 Cx Cy 0;
        0 0 0 -Cy Cx 0;
        0 0 0 0 0 1]
    
    return R
end

function R2d_frame!(element::Element)
    element.R = R2d_frame(element)
end

"""
3d frame rotation matrix
"""
function R3d_frame(element::Element; tol = 1e-4)
    xvec = normalize(element.posEnd .- element.posStart) # local x vector
    CXx, CYx, CZx = xvec # local x cosines

    if norm(cross(xvec, globalY)) < tol #special case for horizontal members aligned with global Y
        Λ = [0. CYx 0.;
            -CYx*cos(element.Ψ) 0 sin(element.Ψ);
            CYx*sin(element.Ψ) 0 cos(element.Ψ)]
    else # all other
        b1 = (-CXx * CYx * cos(element.Ψ) - CZx * sin(element.Ψ)) / sqrt(CXx^2 + CZx^2)
        b2 = sqrt(CXx^2 + CZx^2) * cos(element.Ψ)
        b3 = (-CYx * CZx * cos(element.Ψ) + CXx * sin(element.Ψ)) / sqrt(CXx^2 + CZx^2)

        c1 = (CXx * CYx * sin(element.Ψ) - CZx * cos(element.Ψ)) / sqrt(CXx^2 + CZx^2)
        c2 = -sqrt(CXx^2 + CZx^2) * sin(element.Ψ)
        c3 = (CYx * CZx * sin(element.Ψ) + CXx * cos(element.Ψ)) / sqrt(CXx^2 + CZx^2)

        Λ = [CXx CYx CZx; 
            b1 b2 b3; 
            c1 c2 c3]
    end
    
    R = [Λ zeros(3,9); zeros(3,3) Λ zeros(3,6); zeros(3,6) Λ zeros(3,3); zeros(3,9) Λ]

    return R
end

function R3d_frame!(element::Element)
    element.R = R3d_frame(element)
end

"""
Matches element condition to appropriate rotation matrix
"""
Rdict = Dict(
    (2, :truss) => R2d_truss!,
    (3, :truss) => R3d_truss!,
    (2, :frame) => R2d_frame!,
    (3, :frame) => R3d_frame!
)