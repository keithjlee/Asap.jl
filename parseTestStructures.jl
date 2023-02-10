using Asap, JSON, SparseArrays, LinearAlgebra

f = JSON.parsefile("exampleStructures/biosphere.json")

E = 200e6
G = 77e6
A = 0.001140
Ix = 1.03e6 / 1000^4
J = 2050e3 / 1000^4

sec = Section(A, E, G, Ix, Ix, J)

supportIDs = Vector{Int}()

for sup in f["Supports"]["\$values"]
    push!(supportIDs, sup["NodeInd"] .+ 1)
end

nNodes = length(f["Nodes"]["\$values"])
nodes = Vector{Node}()


for i = 1:nNodes
    node = f["Nodes"]["\$values"][i]
    x = node["Pos"]["X"]
    y = node["Pos"]["Y"]
    z = node["Pos"]["Z"]

    if any(i .== supportIDs)
        node = Node([x, y, z], :fixed)
        node.id = :support
    else
        node = Node([x, y, z], :free)
        node.id = :free
    end

    push!(nodes, node)
end

nElements = length(f["Elems"]["\$values"])
elements = Vector{Element}()

for i = 1:nElements
    elem = f["Elems"]["\$values"][i]
    ids = Int.(elem["NodeInds"]["\$values"] .+ 1)
    push!(elements, Element(nodes, ids, sec))
end

loads = Vector{Load}()
Pvec = [0., 0., -2.0]

for node in nodes[:free]
    load = NodeForce(node, Pvec)
    push!(loads, load)
end

model = Model(nodes, elements, loads)
@time solve!(model; reprocess = true)

#######

f = JSON.parsefile("exampleStructures/DJMM_bridge.json")

ns = f["node_list"]

# elements
es = f["element_list"]

# material properties
mp = f["material_properties"]
E = Float64(mp["youngs_modulus"]) # kn/cm2
A = mp["cross_sec_area"] #cm^2
G = Float64(mp["shear_modulus"]) #kN/cm^2
Iz = mp["Iz"] #cm^4
Iy = mp["Iy"] # cm^4
J = mp["Jx"] #cm^4

sec = Section(A, E, G, Iz, Iy, J)

# Create nodes
NODES = Vector{Node}()

# for DJMM Bridge, flip structure
heights = [n["point"]["Z"] for n in ns]
normalizer = maximum(heights)
fixidx = findall(heights .>= normalizer - 0.001)

for i in eachindex(ns)
    n = ns[i]
    position = [n["point"]["X"], n["point"]["Y"], -n["point"]["Z"] + normalizer] ./ 10
    if any(fixidx .== i)
        fixtype = :fixed
    else
        fixtype = :free
    end
    push!(NODES, Node(position, fixtype))
end

# Create elements
ELEMENTS = Vector{Element}()
for e in es
    indices = e["end_node_ids"] .+ 1
    push!(ELEMENTS, Element(NODES, indices, sec))
    # push!(ELEMENTS, Element(NODES, indices, E, A))
end

# Create self-weight loads
lengths = [e.length for e in ELEMENTS]

# load adjusted to give reasonable analysis/internalforces
Pvec = [0, 0, -0.002]
# Pvec = [0, 0, -0.002]

LOADS = Vector{Load}()
for i in eachindex(NODES)
    if any(i .== fixidx)
        continue
    else
        push!(LOADS, NodeForce(NODES[i],  Pvec))
    end
end

MODEL = Model(NODES, ELEMENTS, LOADS)
@time solve!(MODEL; reprocess = true);