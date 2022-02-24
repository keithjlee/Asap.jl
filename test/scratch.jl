# e1 = [0.0, 0.0, 0.0]
# e2 = [1., 1., 1.]

# CXx, CYx, CZx = e2 .- e1 ./ norm(e2 .- e1)

# Ψ = 0.

# b1 = (-CXx * CYx * cos(Ψ) - CZx * sin(Ψ)) / sqrt(CXx^2 + CZx^2)
# b2 = sqrt(CXx^2 + CZx^2) * cos(Ψ)
# b3 = (-CYx * CZx * cos(Ψ) + CXx * sin(Ψ)) / sqrt(CXx^2 + CZx^2)

# c1 = (CXx * CYx * sin(Ψ) - CZx * cos(Ψ)) / sqrt(CXx^2 + CZx^2)
# c2 = -sqrt(CXx^2 + CZx^2) * sin(Ψ)
# c3 = (CYx * CZx * sin(Ψ) + CXx * cos(Ψ)) / sqrt(CXx^2 + CZx^2)

# Λ = [CXx CYx CZx; 
#     b1 b2 b3; 
#     c1 c2 c3]

# #
# D = sqrt(CXx^2 + CYx^2)
# CXy = -CYx / D
# CYy = CXx / D
# CZy = 0
# CXz = -CXx * CZx / D
# CYz = -CYx * CZx / D
# CZz = D

# Λ2 = [CXx CYx CZx; CXy CYy CZy; CXz CYz CZz]

using Asap
t = 4
# Nodes
a = Node([0., 0., 0.], :frame, :fixed)
b = Node([1., 0., 0.], :frame, :free)
c = Node([1., 1., 0.] / sqrt(2), :frame, :free)
nodes = [a, b]
nodes2 = [a, c]
# elements
begin
    E = 200e6
    A = 9290/1e6
    G = 77e6
    I1 = 113e-6
    I2 = 38.8e-6
    J = 575e-3
end

# Iz is strong axis, Iy is weak axis
e = Element(nodes, [1, 2], E, A, G, I1, I2, J)
elements = [e]

# load
l = Load(nodes, b.position, [0., 0., -10., 0., 0., 0.])
loads = [l]

#structure
s = Structure(nodes, elements, loads)
analyze(s)

u1 = s.U

# Iz is weak axis, Iy is strong axis
e2 = Element(nodes, [1, 2], E, A, G, I2, I1, J)
elements2 = [e2]

s2 = Structure(nodes, elements2, loads)
analyze(s2)

u2 = s2.U

e3 = Element(nodes2, [1,2], E, A, G, I2, I1, J)
elements3 = [e3]
s3 = Structure(nodes2, elements3, loads)
analyze(s3)
u3 = s3.U

Δ_theoretical = -10 * e.length^3 / 3 / E / I1
M_theoretical = 10 * e.length



begin
    element = e
    k1 = element.E * element.A / element.length
    k2 = 12 * element.E * element.Iz / element.length^3
    k3 = 6 * element.E * element.Iz / element.length^2
    k4 = 4 * element.E * element.Iz / element.length
    k5 = 2 * element.E * element.Iz / element.length
    k6 = 12 * element.E * element.Iy / element.length^3
    k7 = 6 * element.E * element.Iy / element.length^2
    k8 = 4 * element.E * element.Iy / element.length
    k9 = 2 * element.E * element.Iy / element.length
    k10 = element.G * element.J / element.length

    #stiffness blocks
    a = [k1 0 0; 0 k2 0; 0 0 k6]
    b = [0 0 0; 0 0 k3; 0 -k7 0]
    c = [k10 0 0; 0 k8 0; 0 0 k4]
    d = [-k10 0 0; 0 k9 0; 0 0 k5]

    k1 = [a b -a b; b' c b d; -a' b' a -b; b' d' -b' c]

    K1 = element.R' * k1 * element.R

    E = element.E
    length = element.length
    A = element.A
    Iz = element.Iz
    Iy = element.Iy
    G = element.G
    J = element.J
    k2 = E / length^3 * [
        A*length^2 0 0 0 0 0 -A*length^2 0 0 0 0 0;
        0 12Iz 0 0 0 6length*Iz 0 -12Iz 0 0 0 6length*Iz;
        0 0 12Iy 0 -6length*Iy 0 0 0 -12Iy 0 -6length*Iy 0;
        0 0 0 G*J*length^2/E 0 0 0 0 0 -G*J*length^2/E 0 0;
        0 0 -6length*Iy 0 4length^2*Iy 0 0 0 6length*Iy 0 2length^2*Iy 0;
        0 6length*Iz 0 0 0 4length^2*Iz 0 -6length*Iz 0 0 0 2length^2*Iz;
        -A*length^2 0 0 0 0 0 A*length^2 0 0 0 0 0;
        0 -12Iz 0 0 0 -6length*Iz 0 12Iz 0 0 0 -6length*Iz;
        0 0 -12Iy 0 6length*Iy 0 0 0 12Iy 0 6length*Iy 0;
        0 0 0 -G*J*length^2/E 0 0 0 0 0 G*J*length^2/E 0 0;
        0 0 -6length*Iy 0 2length^2*Iy 0 0 0 6length*Iy 0 4length^2*Iy 0;
        0 6length*Iz 0 0 0 2length^2*Iz 0 -6length*Iz 0 0 0 4length^2*Iz
        ]

    K2 = element.R' * k2 * element.R
end