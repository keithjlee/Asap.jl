# 3 main data structures: Node, Element, Structure

# defining all nodes: Node([position], type, boundary conditions)
node1 = Node([72.0, 0.0, 0.0], :truss, :yfixed)
node2 = Node([0.0, 36.0, 0.0], :truss, :fixed)
node3 = Node([0.0, 36.0, 72.0], :truss, :fixed)
node4 = Node([0.0, 0.0, -48.0], :truss, :fixed)

# The following is equivalent: Node([position], [DOFs])
node1 = Node([72.0, 0.0, 0.0], [true, false, true])
node2 = Node([0.0, 36.0, 0.0], [false, false, false])
node3 = Node([0.0, 36.0, 72.0], [false, false, false])
node4 = Node([0.0, 0.0, -48.0], [false, false, false])

# group nodes together
nodes = [node1, node2, node3, node4]

# Material properties
E = 1.2e6

# defining elements
e1 = Element(nodes, [1,2], E, 0.302)
e2 = Element(nodes, [1,3], E, 0.729)
e3 = Element(nodes, [1,4], E, 0.187)

# group elements
elements = [e1, e2, e3]

# loads
l1 = Load(nodes, node1.position, [0.0,0.0,-1000.0])

# group loads (in this case just 1)
loads = [l1]

# assemble structure
structure = Structure(nodes, elements, loads)

# analyze structure
analyze(structure)

# relevant information is now stored in the Structure datatype, e.g.
displacements = structure.U

# the Geometry type stores things like displaced positions for visualization
geo = Geometry(structure)

# structurePlot is used to plot the Geometry type
fig = structurePlot(geo)