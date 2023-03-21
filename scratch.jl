using Asap, JSON, SparseArrays, LinearAlgebra, Statistics, kjlMakie, GLMakie


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

#### Simply Supported Beam


#cantilever


dpos = Asap.disp(model, e, 40, 100)

lines(Point3.(dpos), axis = (aspect = DataAspect(),))

#cantilever middle load
P=1000.; 
m = 200.;
L=144.; 
E=30e6; 
Is=57.; 

a = .625;

n1 = Node([0., 0., 0.], :fixed)
n2 = Node([L, 0., 0.], :yfixed)

nodes = [n1, n2]

sec = Section(1., E, 1., Is, Is, 1.)

e = Element(nodes, [1,2], sec)
elements = [e]

p1 = LineLoad(e, [0., -m, 0.])
p2 = PointLoad(e, a, [0., -P, 0.])

loads = [p1, p2]

model = Model(nodes, elements, loads)
planarize!(model)
solve!(model; reprocess = true)

dpos = Asap.disp(model, e, 40, 100)

p1 = Point3.(dpos)

dan(x) = 3.7229e-7 * x^3 - 5.361e-5 * x^2

xr = range(0, e.length, 40)

p2 = Point3.([[x, dan(x) * 100, 0.] for x in xr])

begin
    fig = Figure()
    ax = Axis(fig[1,1])

    lines!(p1, color = blue)

    lines!(p2, color = :black)
    fig
end

b = 1 - a

P * a^2 * b / L^2 + m * L^2 / 12