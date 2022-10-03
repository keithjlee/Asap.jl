# 3D frame test: Example 8.3 Ferreira "MATLAB Codes for Finite Element Analysis"
# in kN, m
using Asap
# nodes
n1 = Node([0., 0., 0.], :frame, :free)
n2 = Node([3., 0., 0.], :frame, :fixed)
n3 = Node([0., 0., -3.], :frame, :fixed)
n4 = Node([0., -4., 0.], :frame, :fixed)

nodes = [n1, n2, n3, n4]

# elements
E = 210e6
A = 0.02
Iy = 10e-5
Iz = 20e-5
J = 5e-5
G = 84e6

e1 = Element(nodes, [1,2], E, A, G, Iz, Iy, J)
e2 = Element(nodes, [1,3], E, A, G, Iz, Iy, J)
e3 = Element(nodes, [1,4], E, A, G, Iz, Iy, J)

elements = [e1, e2, e3]

# loads
l = Load(nodes, n1.position, [-10., 0., 20., 0., 0., 0.])

loads = [l]

ex83 = Structure(nodes, elements, loads)
analyze!(ex83)

d_textbook = [-7.05e-6, -7e-8, 1.418e-5, 1.45e-6, 1.75e-6, 1.14e-6]

@test round.(d_textbook, digits = 8) ≈ round.(n1.disp, digits = 8)
