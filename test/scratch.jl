using Asap

# in kN, m

# nodes
n1 = Node([0., 0.], :frame, :fixed)
n2 = Node([0., 10.], :frame, :free)
n3 = Node([8., 10.], :frame, :pinfixed)

nodes = [n1, n2, n3]

# elements
E = 200e6 #kN/m^2
A = 0.00474 #m^2
I = 0.0000222 #m^4

e1 = Element(nodes, [1,2], E, A, I)
e2 = Element(nodes, [2,3], E, A, I)

elements = [e1, e2]

# loads
l1 = Load(nodes, n2.position, [120., -37.5, 125.])
l2 = Load(nodes, n3.position, [0., 0., 75.])

loads = [l1, l2]

ex65 = Structure(nodes, elements, loads)
analyze(ex65)
ex65.U