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

testDict = Dict(
    (2, :truss, :free) => [true, true],
    (2, :truss, :fixed) => [false, false],
    (2, :truss, :xfree) => [true, false],
    (2, :truss, :yfree) => [false, true],
    (2, :frame, :free) => [true, true, true],
    (2, :frame, :pinfixed) => [false, false, true],
    (2, :frame, :fixed) => [false, false, false],
    (2, :frame, :xfree) => [true, false, true],
    (2, :frame, :yfree) => [false, true, false],
    (3, :truss, :free) => [true, true, true],
    (3, :truss, :fixed) => [false, false, false],
    (3, :truss, :zfixed) => [true, true, false],
    (3, :truss, :xfixed) => [false, true, true],
    (3, :truss, :yfixed) => [true, false, true],
    (3, :truss, :xfree) => [true, false, false],
    (3, :truss, :yfree) => [false, true, false],
    (3, :truss, :zfree) => [false, false, true],
    (3, :frame, :free) => [true, true, true, true, true, true],
    (3, :frame, :pinfixed) => [false, false, false, true, true, true],
    (3, :frame, :fixed) => [false, false, false, false, false, false],
    (3, :frame, :xfixed) => [false, true, true, true, true, true],
    (3, :frame, :yfixed) => [true, false, true, true, true, true],
    (3, :frame, :zfixed) => [true, true, false, true, true, true]
)