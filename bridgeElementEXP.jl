using LinearAlgebra, kjlMakie, Asap, SteelSections
set_theme!(kjl_light)

w = W("W310X97")
# w = HSSRound("HSS88.9X3.2")
mat = Steel_Nmm
sec = toASAPframe(w, mat.E, mat.G)
secBig = toASAPframe(W("W460X97"), mat.E, mat.G)

# definitions
begin
    n1 = Node([0., 0., 0.], :pinned)
    n2 = Node([10_000., 0., 0.], :zfixed)
    n3 = Node([0., 4_000., 0.], :pinned)
    n4 = Node([10_000., 4_000., 0.], :zfixed)

    nodes = [n1, n2, n3, n4]

    e1 = Element(n1, n2, secBig)
    e2 = Element(n3, n4, secBig)
    e3 = Element(n1, n4, secBig)

    girders = [e1, e2]

    n = 20

    joists = [BridgeElement(e1, x, e2, x, sec) for x in range(0, 1, n)[2:end-1]]

    # b1 = Asap.BridgeElement(e1, 0.5, e2, 0.5, sec)
    # b2 = Asap.BridgeElement(e1, 0.2, e2, 0.2, sec)
    # b3 = Asap.BridgeElement(e2, 0.6, e1, 0.1, sec)


    # elements = [[e1, e2, e3] ; [b1, b2, b3]]

    elements = [girders; joists]

    gravity = [LineLoad(j, [0., 0., -10.]) for j in joists]

    uplifts = [LineLoad(rand(joists), [0., 0., 40.]) for _ = 1:4]


    loads = [gravity; uplifts]

    # l1 = LineLoad(e2, [0., 0., -2.])
    # l2 = PointLoad(b1, 0.5, [0., 0., -10])
    # l3 = PointLoad(e1, 0.3, [0., 0., -3.])
    # l4 = LineLoad(e3, [0., 0., -1.])
    # l5 = LineLoad(b3, [0., 0., 4.])

    # loads = [l1, l2, l3, l4, l5]


end

## walkthrough process_bridge!
model = Model(nodes, elements, loads)

@time Asap.process_bridge!(model);

solve!(model)

nu = Point3.([n.position for n in model.nodes])
eu = vcat([nu[e.nodeIDs] for e in model.elements]...)

fac = 50
n = 10

nd = Point3.([n.position .+ fac * n.displacement[1:3] for n in model.nodes])
ed = [Point3.(eachcol(displacedshape(e; factor = fac, n = n))) for e in model.elements]

begin
    fig = Figure()
    ax = Axis3(fig[1,1],
        # aspect = (1,1,1),
        aspect = :data
        )

    e1 = linesegments!(eu,
        color = (:black, 0.5))

    n1 = scatter!(nu,
        color = :black)

    e2 = lines!.(ed,
        color = :black,
        linewidth = 3)

    n2 = scatter!(nd,
        color = :black)

    fig
end