"""
Transformation matrix of an element
"""
function R(element::Element; tol = 1e-4)
    
    xvec = local_x(element)
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
Transformation matrix of a directional vector and roll angle
"""
function R(xvec::Vector{Float64}, Ψ; tol = 1e-4)

    CXx, CYx, CZx = normalize(xvec) # local x cosines

    cΨ = cos(Ψ)
    sΨ = sin(Ψ)


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
Truss transformation matrix
"""
function R(element::TrussElement)
    Cx, Cy, Cz = local_x(element)
    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]
end

function R(vec::Vector{Float64})
    Cx, Cy, Cz = normalize(vec)
    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]
end

function R(vec::SubArray)
    Cx, Cy, Cz = normalize(vec)
    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]
end

