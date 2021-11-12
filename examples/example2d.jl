using Asap

# Generate nodes
structureType = :truss
begin
    a = Node([0.0, 0.0], structureType, :fixed)
    b = Node([7.5, 10.0], structureType, :free)
    c = Node([15.0, 10.0], structureType, :free)
    d = Node([22.5, 10.0], structureType, :free)
    e = Node([30.0, 0.0], structureType, :xfree)
    f = Node([22.5, 0.0], structureType, :free)
    g = Node([15.0, 0.0], structureType, :free)
    h = Node([7.5, 0.0], structureType, :free)
    nodes = [a,b,c,d,e,f,g,h]
end

# HSS141x6.4 Round sections material properties
begin
    E = 200e6 #kN/m^2
    A = 2440 / 1e6 #m^2
    G = 77e6 #kN/m^2
    Iz = 5.61e-6 #m^4
    Iy = 5.61e-6 #m^4
    J = 11.2e-6 #m^4
end

# generate elements
if structureType == :truss
    matProps = (E,A)
else
    matProps = (E,A,Iz)
end

begin
    e1 = Element(nodes, [1,2], matProps...)
    e2 = Element(nodes, [2,3], matProps...)
    e3 = Element(nodes, [3,4], matProps...)
    e4 = Element(nodes, [4,5], matProps...)
    e5 = Element(nodes, [5,6], matProps...)
    e6 = Element(nodes, [6,7], matProps...)
    e7 = Element(nodes, [7,8], matProps...)
    e8 = Element(nodes, [8, 1], matProps...)
    e9 = Element(nodes, [2, 8], matProps...)
    e10 = Element(nodes, [2,7], matProps...)
    e11 = Element(nodes, [3,7], matProps...)
    e12 = Element(nodes, [4,7], matProps...)
    e13 = Element(nodes, [4,6], matProps...)
    elements = [e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13]
end

# create loads
if structureType == :truss
    forces = [0.0, -20.0]
else
    forces = [0.0, -20.0, 0.0]
end

begin
    p1 = Load(nodes, nodes[6].position, forces)
    p2 = Load(nodes, nodes[7].position, forces)
    p3 = Load(nodes, nodes[8].position, forces)

    loads = [p1, p2, p3]
end

# assemble structure
truss = Structure(nodes, elements, loads)

# analyze
@time analyze(truss)

# output geometry
geo = Geometry(truss)

# display figure (this is quite crude at the moment)
# scaleToForce scales line thicknesses to axial force values
fig = structurePlot(geo; scaleToForce = true)


# 3 main data structures: Node, Element, Structure
# Look in types.jl for the data fields and different function methods of creating each type

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