using Asap, JSON, SparseArrays, LinearAlgebra, Statistics, GLMakie, kjlMakie, Zygote


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
structure = TrussModel(nodes, elements, loads);
@time solve!(structure; reprocess = true)

################
# NAIVE SHAPE OPTIMIZATION
################
###CONSTANTS
#elemental stiffness matrices in LCS
Kstore = Asap.localK.(structure.elements)
#free degrees of freedom
freedof = structure.freeDOFs
#load vector
P = structure.P
#Connectivity Matrix
C = Asap.connectivity(structure)
#Initial positions
X = Asap.nodePositions(structure)
#Row/Column vectors for global K assembly
I = Vector{Int64}()
J = Vector{Int64}()

for element in structure.elements
    idx = element.globalID
    @inbounds for j = 1:6
        @inbounds for i = 1:6
            push!(I, idx[i])
            push!(J, idx[j])
        end
    end
end

@time begin
    x0 = X[:, 3]

    Xvec = [X[:, 1:2] x0]

    E = normalize.(eachrow(C * Xvec))

    Rs = Asap.R.(E)

    KeGlobal = transpose.(Rs) .* Kstore .* Rs;
    V = vcat(vec.(KeGlobal)...)

    # V = vcat([vec(r' * k * r) for (r,k) in zip(Rs, Kstore)]...)

    K = sparse(I, J, V)

    U = K[freedof, freedof] \ P[freedof]

    O = U' * P[freedof]
end

Cd = Matrix(C)

Is = Vector{Vector{Int64}}()
Js = Vector{Vector{Int64}}()

for element in structure.elements
    idx = element.globalID
    ii = Vector{Int64}()
    jj = Vector{Int64}()

    for j = 1:6
        for i = 1:6
            push!(ii, idx[i])
            push!(jj, idx[j])
        end
    end

    push!(Is, ii)
    push!(Js, jj)
end

n = structure.nDOFs

ks = [sparse(i, j, vec(k), n, n) for (i, j, k) in zip(Is, Js, KeGlobal)]

function obj(x)
    Xvec = [X[:, 1:2] x]

    E = normalize.(eachrow(C * Xvec))

    Rs = Asap.R.(E)

    #OPTION 1: FASTER BUT UNSTABLE
    KeGlobal = vec.(transpose.(Rs) .* Kstore .* Rs)
    # V = vcat(KeGlobal...)
    # K = sparse(I, J, V)

    #OPTION 2: MEMORY INTENSIVE BUT WORKS
    ks = [sparse(i, j, k, n, n) for (i, j, k) in zip(Is, Js, KeGlobal)]
    K = sum(ks)

    U = K[freedof, freedof] \ P[freedof]

    O = U' * P[freedof]

    return O
end


x0 = X[:, 3] .+ rand(structure.nNodes)

@time obj(x0)

@time g = gradient(x-> obj(x), x0)[1];

