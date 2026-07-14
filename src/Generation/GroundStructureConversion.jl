function to_truss(gs::GroundStructure, section::AbstractSection; load = [1., 0., 0.], support = :xy, planar = true)

    @assert in(support, [:xy, :x, :y, :corner])

    nodes = [Node([xy..., 0.], :free, :free) for xy in gs.xy]

    if support == :xy
        isupport = [gs.igrid[[1,end],:][:]; gs.igrid[2:end-1, [1,end]][:]]
    elseif support == :x
        isupport = gs.igrid[[1,end], :][:]
    elseif support == :y
        isupport = gs.igrid[:, [1,end]][:]
    else
        isupport = [gs.igrid[1,1], gs.igrid[1,end], gs.igrid[end,1], gs.igrid[end,end]]
    end

    for node in nodes[isupport]
        fixnode!(node, :pinned)
        node.id = :support
    end 

    elements = [TrussElement(nodes[id]..., section, :free) for id in gs.ielements]

    loads = [NodeForce(node, load) for node in nodes[:free]]

    model = Model(nodes, elements, loads)

    planar && (planarize!(model))

    solve!(model)

    return model

end

function to_frame(gs::GroundStructure, section::Section; load = [0., 0., -1.], support = :xy, planar = false, support_type = :pinned)

    @assert in(support, [:xy, :x, :y, :corner])
    @assert in(support_type, [:pinned, :fixed])

    nodes = [Node([xy..., 0.], :free, :free) for xy in gs.xy]

    if support == :xy
        isupport = [gs.igrid[[1,end],:][:]; gs.igrid[2:end-1, [1,end]][:]]
    elseif support == :x
        isupport = gs.igrid[[1,end], :][:]
    elseif support == :y
        isupport = gs.igrid[:, [1,end]][:]
    else
        isupport = [gs.igrid[1,1], gs.igrid[1,end], gs.igrid[end,1], gs.igrid[end,end]]
    end

    for node in nodes[isupport]
        fixnode!(node, support_type)
        node.id = :support
    end 

    elements = [FrameElement(nodes[id]..., section, :free) for id in gs.ielements]

    loads = [NodeForce(node, load) for node in nodes[:free]]

    model = Model(nodes, elements, loads)

    planar && (planarize!(model))

    solve!(model)

    return model

end