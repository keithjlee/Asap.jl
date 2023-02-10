using Asap, JSON, SparseArrays, LinearAlgebra, Statistics, GLMakie, kjlMakie


E = 200e6 #kN/m^2
density = 80 #kN/m^3

filename = "exampleStructures/space_truss_00023.json"
truss = JSON.parsefile(filename)

# nodes and elements
ns = truss["node_list"]
es = truss["element_list"]

# Create nodes
nodes = Vector{TrussNode}()
for n in ns
    position = [n["point"]["X"], n["point"]["Y"], n["point"]["Z"]]
    if n["is_grounded"] == 1
        node = TrussNode(position, :fixed)
        node.id = :support
    else
        node = TrussNode(position, :free)
        node.id = :free
    end

    push!(nodes, node)
end

# Create elements
elements = Vector{TrussElement}()

sec = TrussSection(0.002, E)
for e in es
    indices = Int.(e["end_node_ids"] .+ 1)
    push!(elements, TrussElement(nodes, indices, sec))
end

# Create self-weight loads
lengths = [e.length for e in elements]
Pvec = [0, 0, -3.0] # 1kN downwards force 
loads = Vector{NodeForce}()
for node in nodes[:free]
    push!(loads, NodeForce(node, Pvec))
end
structure = TrussModel(nodes, elements, loads)
@time solve!(structure; reprocess = true)