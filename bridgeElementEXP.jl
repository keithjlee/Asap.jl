using LinearAlgebra, kjlMakie, Asap, SteelSections, Random
set_theme!(kjl_light)

w = W("W310X97")
# w = HSSRound("HSS88.9X3.2")
mat = Steel_Nmm
sec = toASAPframe(w, mat.E, mat.G)
secBig = toASAPframe(W("W460X97"), mat.E, mat.G)

# definitions
begin
    n1 = Node([0., 0., 0.], :pinned)
    n2 = Node([10_000., 0., 0.], :xfree)
    n3 = Node([0., 4_000., 0.], :pinned)
    n4 = Node([10_000., 4_000., 0.], :xfree)
    n5 = Node([5_000., 5_000., 1300.], :pinned)

    nodes = [n1, n2, n3, n4, n5]

    e1 = Element(n1, n2, secBig)
    e2 = Element(n3, n4, secBig)
    e3 = Element(n3, n5, secBig)
    e4 = Element(n4, n5, secBig)

    girders = [e1, e2, e3, e4]
    for e in girders 
        e.id = :girders
    end

    njoists = 100
    joists1 = [BridgeElement(e1, x, e2, x, sec) for x in range(0, 1, njoists)[2:end-1]]

    for e in joists1
        e.id = :joists_horizontal
    end

    njoists2 = 50
    joists2 = [BridgeElement(e2, x/2, e3, x, sec) for x in range(0, 1, njoists2)[2:end-1]]
    joists3 = [BridgeElement(e2, 1 - x/2, e4, x, sec) for x in range(0, 1, njoists2)[2:end-1]]


    elements = [girders; joists1; joists2; joists3]

    gravity = [LineLoad(j, [0., 0., -0.5]) for j in joists1]

    upwardForces = [PointLoad(j, 0.5, [0., 0., 75.]) for j in [joists2; joists3]]


    loads = [gravity; upwardForces]

    push!(loads, PointLoad(e1, 0.5, [0., 0., 5000.]))

    # l1 = LineLoad(e2, [0., 0., -2.])
    # l2 = PointLoad(b1, 0.5, [0., 0., -10])
    # l3 = PointLoad(e1, 0.3, [0., 0., -3.])
    # l4 = LineLoad(e3, [0., 0., -1.])
    # l5 = LineLoad(b3, [0., 0., 4.])

    # loads = [l1, l2, l3, l4, l5]}

    shuffle!(elements)


end;

## walkthrough process_bridge!
model = Model(nodes, elements, loads);
@time Asap.process_bridge!(model);


solve!(model)
begin
    nu = Point3.([n.position for n in model.nodes])
    eu = vcat([nu[e.nodeIDs] for e in model.elements]...)

    fac = Observable(50)
    n = 10

    # nd = Point3.([n.position .+ fac * n.displacement[1:3] for n in model.nodes])
    ed = [@lift(Point3.(eachcol(displacedshape(e; factor = $fac, n = n)))) for e in model.elements]
    # edc = [repeat([e.elementID], length(d)) for (e,d) in zip(model.elements, ed)]
end

begin
    fig = Figure()
    ax = Axis3(fig[1,1],
        # aspect = (1,1,1),
        aspect = :data
        )

    e1 = linesegments!(eu,
        color = blue)

    # n1 = scatter!(nu,
    #     color = :black)

    e2 = lines!.(ed,
        color = :black,
        # color = edc,
        # colorrange = (1,20),
        # colormap = :tempo,
        linewidth = 3)

    on(fac) do _
        reset_limits!(ax)
    end

    # n2 = scatter!(nd,
    #     color = :black)

    fig
end


