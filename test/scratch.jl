using Asap
n1 = Node([0., 0.] .* 12, :frame, :fixed)
n2 = Node([10., 20.] .* 12, :frame, :free)
n3 = Node([30., 20.] .* 12, :frame, :fixed)

nodes = [n1, n2, n3]

# elements
E = 29e3 #ksi
A = 11.8 # in^2
I = 310. # in^4

e1 = Element(nodes, [1,2], E, A, I)
e2 = Element(nodes, [2,3], E, A, I)

elements = [e1, e2]

# loads
l1 = Load(nodes, n2.position, [0., -60., -750.])

loads = [l1]

# assembly + analysis
ex66 = Structure(nodes, elements, loads)
analyze!(ex66)

f1 = Rdict[(2, :truss)]
