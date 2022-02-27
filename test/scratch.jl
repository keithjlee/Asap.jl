using Asap
# Nodes
a = Node([0., 0., 0.], :frame, :fixed)
b = Node([0., 1., 0.], :frame, :free)
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
r1 = s.reactions

# Iz is weak axis, Iy is strong axis
e2 = Element(nodes, [1, 2], E, A, G, I2, I1, J)
e2.Ψ = 0.
elements2 = [e2]

s2 = Structure(nodes, elements2, loads)
analyze(s2)

u2 = s2.U
r2 = s2.reactions