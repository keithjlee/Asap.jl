e1 = [0.0, 0.0, 0.0]
e2 = [1., 1., 1.]

CXx, CYx, CZx = e2 .- e1 ./ norm(e2 .- e1)

Ψ = 0.

b1 = (-CXx * CYx * cos(Ψ) - CZx * sin(Ψ)) / sqrt(CXx^2 + CZx^2)
b2 = sqrt(CXx^2 + CZx^2) * cos(Ψ)
b3 = (-CYx * CZx * cos(Ψ) + CXx * sin(Ψ)) / sqrt(CXx^2 + CZx^2)

c1 = (CXx * CYx * sin(Ψ) - CZx * cos(Ψ)) / sqrt(CXx^2 + CZx^2)
c2 = -sqrt(CXx^2 + CZx^2) * sin(Ψ)
c3 = (CYx * CZx * sin(Ψ) + CXx * cos(Ψ)) / sqrt(CXx^2 + CZx^2)

Λ = [CXx CYx CZx; 
    b1 b2 b3; 
    c1 c2 c3]

#
D = sqrt(CXx^2 + CYx^2)
CXy = -CYx / D
CYy = CXx / D
CZy = 0
CXz = -CXx * CZx / D
CYz = -CYx * CZx / D
CZz = D

Λ2 = [CXx CYx CZx; CXy CYy CZy; CXz CYz CZz]