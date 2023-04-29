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

    e1 = Element(n1, n2, sec)
    e2 = Element(n3, n4, sec)
    e3 = Element(n1, n4, sec)

    b1 = Asap.BridgeElement(e1, 0.5, e2, 0.5, sec)
    b2 = Asap.BridgeElement(e1, 0.2, e2, 0.2, sec)
    b3 = Asap.BridgeElement(e2, 0.6, e1, 0.1, sec)

    bes = [Asap.BridgeElement(e1, rand(), e2, rand(), sec) for _ = 1:20]

    elements = [[e1, e2, e3] ; bes]

    l1 = LineLoad(e2, [0., 0., -3])
    l2 = PointLoad(rand(bes), 0.5, [0., 0., -10e3])
    l3 = PointLoad(e1, 0.3, [0., 0., 3.])
    l4 = LineLoad(e3, [0., 0., -1.])

    loads = [l1, l2, l3, l4]
end

#block 1: identifiers
Asap.makeids!(nodes)
Asap.makeids!(elements)
Asap.makeids!(loads)

# nodes
begin
    n_dof = 6
    dofset = collect(0:n_dof-  1)
    for (i, node) in enumerate(nodes)
        node.globalID = i * n_dof - (n_dof - 1) .+ dofset
    end
end

# loads
begin
    for load in loads
        Asap.assign!(load)
    end

end

# elements
newPos = Vector{Vector{Float64}}()

for element in elements
    if typeof(element) == Element
        element.Q = zeros(12) # reset Qf
        element.R = Asap.R(element)
        element.LCS = Asap.lcs(element, element.Ψ)
        element.length = Asap.dist(element.nodeStart, element.nodeEnd)
        Asap.makeK!(element)

        element.nodeIDs = [element.nodeStart.nodeID, element.nodeEnd.nodeID]

        idStart = element.nodeStart.globalID
        idEnd = element.nodeEnd.globalID

        element.globalID = [idStart; idEnd]
    else
        es = element.elementStart
        npos = es.nodeStart.position .+ es.LCS[1] .* es.length .* element.posStart

        ee = element.elementEnd
        npos2 = ee.nodeStart.position .+ ee.LCS[1] .* ee.length .* element.posEnd

        push!(newPos, npos, npos2)
    end
end

iElements = findall(typeof.(elements) .== Element)
iBridges = findall(typeof.(elements) .== Asap.BridgeElement)

n = length(iElements)
m = length(iBridges)

BMmat = zeros(Int64, n, m)
Nodemat = [getproperty.(elements[iElements], :nodeStart) Matrix{Node}(undef, n, m) getproperty.(elements[iElements], :nodeEnd)]
ordermat = [zeros(n, m+1) ones(n)]
idvec = getproperty.(elements[iElements], :elementID)

i = 1

# startposStore = Vector{Float64}()
# endposStore = Vector{Float64}()

g2lID = Dict(getproperty.(elements[typeof.(elements) .== Element], :elementID) .=> collect(1:n))

for e in elements
    if typeof(e) == Asap.BridgeElement
        # push!(startposStore, e.posStart)
        # push!(endposStore, e.posEnd)

        rowstart = g2lID[e.elementStart.elementID]
        rowend = g2lID[e.elementEnd.elementID]

        BMmat[rowstart, i] = -1
        BMmat[rowend, i] = 1

        ordermat[rowstart, i + 1] = e.posStart
        ordermat[rowend, i + 1] = e.posEnd

        pos1 = e.elementStart.nodeStart.position .+ e.elementStart.LCS[1] * e.posStart * e.elementStart.length

        Nodemat[rowstart, i + 1] = Node(pos1, :free)

        pos2 = e.elementEnd.nodeStart.position .+ e.elementEnd.LCS[1] * e.posEnd * e.elementEnd.length

        Nodemat[rowend, i + 1] = Node(pos2, :free)

        i += 1
    end
end

iActive = findall(sum.(eachrow(BMmat)) .!= 0)
BMmat = BMmat[iActive,:]
Nodemat = Nodemat[iActive]

newEls = Vector{Element}()
#generate shattered elements
for (j, row, order) in zip(1:size(Nodemat)[1], eachrow(Nodemat), eachrow(ordermat))
    ordernodes = row[sortperm(order)]
    shattered =[Element(ordernodes[i], ordernodes[i+1], e1.section) for i = 1:length(row)-1]

    for e in shattered
        e.elementID = j
    end
    push!(newEls, shattered...)
end

for e in newEls
    e.length = Asap.dist(e.nodeStart, e.nodeEnd)
end

#generate bridge elements
for (i, bm) in enumerate(elements[typeof.(elements) .== Asap.BridgeElement])
    order = sortperm(BMmat[:, i])

    el = Element(Nodemat[order, i+1]..., bm.section)
    el.Ψ = bm.Ψ
    el.id = bm.id
    el.elementID = bm.elementID

    push!(newEls, el)
end


#transfer loads
newLoads = Vector{Asap.Load}()

#for distributed loads
eid1 = l1.element.elementID
itransfer = findall(getproperty.(newEls, :elementID) .== eid1)
for i in itransfer
    push!(newLoads, LineLoad(newEls[i], l1.value))
end

# for point load on bridge element
eid2 = l2.element.elementID
itransfer = findfirst(getproperty.(newEls, :elementID) .== eid2)
push!(newLoads, PointLoad(newEls[itransfer], l2.position, l2.value))

#for point load on an initial beam
eid3 = l3.element.elementID
loadpos = l3.element.length * l3.position
itransfer = findall(getproperty.(newEls, :elementID) .== eid3)
endpositions = cumsum(getproperty.(newEls[itransfer], :length))
ielement = findfirst(endpositions .> loadpos)
newfrac = (loadpos - endpositions[ielement - 1]) / (endpositions[ielement] - endpositions[ielement - 1])

push!(newLoads, PointLoad(newEls[itransfer[ielement]], newfrac, l3.value))



newNodes = vec(Nodemat)

newModel = Model(newNodes, newEls, newLoads)
solve!(newModel; reprocess = true)

pos = Point3.([n.position for n in newModel.nodes])
undisp = vcat([pos[e.nodeIDs] for e in newModel.elements]...)

fac = 1
disp = [Point3.(eachcol(displacedshape(e, factor = fac))) for e in newModel.elements]

begin
    fig = Figure()
    ax=  Axis3(fig[1,1],
        # aspect = :data,
        aspect = (1,1,1)
        )

    ud = linesegments!(undisp,
        color = :black)

    d = lines!.(disp,
        color = :black,
        linewidth = 3)

    fig
end


begin
    n1 = Node([0., 0., 0.], :pinned)
    n2 = Node([2000., 0., 0.], :free)
    n3 = Node([5000., 0., 0.], :free)
    n4 = Node([10000., 0., 0.], :zfixed)

    offset = [0., 4000., 0.]
    n5 = Node(n1.position .+ offset, :pinned)
    n6 = Node(n2.position .+ offset, :free)
    n7 = Node(n3.position .+ offset, :free)
    n8 = Node(n4.position .+ offset, :zfixed)

    nodes = [n1, n2, n3, n4, n5, n6, n7, n8]

    e1 = Element(n1, n2, secBig)
    e2 = Element(n2, n3, secBig)
    e3 = Element(n3, n4, secBig)

    e4 = Element(n5, n6, secBig)
    e5 = Element(n6, n7, secBig)
    e6 = Element(n7, n8, secBig)

    e7 = Element(n2, n6, sec)
    e8 = Element(n3, n7, sec)

    elements = [e1, e2, e3, e4, e5, e6, e7, e8]

    loads = [LineLoad(e, [0., 0., -2.]) for e in elements]
    model = Model(nodes, elements, loads)
    solve!(model)

end

newloadset = [NodeForce(n, [0., 0., -20e3 * rand()]) for n in model.nodes if all(n.dof)]
solve!(model, newloadset)

pos = Point3.([n.position for n in model.nodes])
undisp = vcat([pos[e.nodeIDs] for e in model.elements]...)

fac = 100
disp = [Point3.(eachcol(displacedshape(e, factor = fac))) for e in model.elements]

begin
    fig = Figure()
    ax=  Axis3(fig[1,1],
        aspect = :data,
        # aspect = (1,1,1)
        )

    ud = linesegments!(undisp,
        color = :black)

    d = lines!.(disp,
        color = :black,
        linewidth = 3)

    fig
end