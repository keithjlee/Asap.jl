using Asap, JSON, Statistics

# Taken from NCCR compas.dev example files
test1 = JSON.parsefile("examples/compas_fea_truss_frame.json")

# nodes (in mm)
ns = test1["node_list"]

# elements
es = test1["element_list"]

# material properties
mp = test1["material_properties"]
E = Float64(mp["youngs_modulus"]) # kn/cm2
A = mp["cross_sec_area"] #cm^2
G = Float64(mp["shear_modulus"]) #kN/cm^2
Iz = mp["Iz"] #cm^4
Iy = mp["Iy"] # cm^4
J = mp["Jx"] #cm^4

# Create nodes
nodes = Vector{Node}()

for i = 1:length(ns)
    n = ns[i]
    position = [n["point"]["X"], n["point"]["Y"], n["point"]["Z"]] ./ 10

    if n["is_grounded"] > 0
        fixtype = :fixed
    else
        fixtype = :free
    end
    push!(nodes, Node(position, :frame, fixtype))
end

# Create elements
elements = Vector{Element}()
for e in es
    indices = e["end_node_ids"] .+ 1
    push!(elements, Element(nodes, indices, E, A, G, Iz, Iy, J))
end

# Create self-weight loads
lengths = [e.length for e in elements]
P = mean(lengths) * A * mp["density"]/100^3 / 2
Pvec = [0, 0, -P, 0, 0, 0]

loads = Vector{Load}()
for node in nodes
    push!(loads, Load(nodes, node.position, Pvec))
end

structure = Structure(nodes, elements, loads)

analyze(structure)

geo = Geometry(structure)
structurePlot(geo; scaleToForce = true)
