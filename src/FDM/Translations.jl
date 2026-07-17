"""
    to_network(model::Model)

Convert a solved (truss) model into an equivalent FDM Network.
"""
function to_network(model::Model)
    if isnothing(model.results)
        error("Analyze truss model before conversion")
    end

    # convert nodes
    nodeset = Vector{FDMnode{Float64}}()

    for node in model.nodes
        pos = Vector(node.position)
        id = node.id

        #any structural support becomes a FULL anchor. Deliberately NOT
        #per-axis: partially-anchored directions with mixed-sign force
        #densities (a converted truss has compression members) are singular
        #force-density systems — e.g. a roller's free direction has no
        #anchor and no load, so its geometry is indeterminate. Build
        #networks with per-axis `FDMnode` fixity directly when you mean it.
        fdmn = FDMnode(pos, all(node.fixity[1:3]))
        fdmn.id = id

        push!(nodeset, fdmn)
    end

    #convert loads
    loadset = Vector{FDMload{Float64}}()

    for load in model.loads
        i = load.node.index

        push!(loadset, FDMload(nodeset, i, Vector(load.value)))
    end

    #convert elements
    elset = Vector{FDMelement{Float64}}()

    for element in model.elements
        istart, iend = element.nodeStart.index, element.nodeEnd.index
        id = element.id
        q = axial_force(model.results, element) / length(element)

        el = FDMelement(nodeset[[istart, iend]]..., q)
        el.id = id

        push!(elset, el)
    end

    network = Network(nodeset, elset, loadset)
    solve!(network)

    return network
end

"""
    to_truss(network::Network, section::AbstractSection)

Convert a solved FDM Network into an equivalent truss model with a given section. All fixed nodes are converted into pinned boundary conditions.
"""
function to_truss(network::Network, section::AbstractSection)
    if network.cache === nothing
        error("Analyze network before conversion")
    end

    nodeset = Vector{Node{Float64}}()
    elset = Vector{AbstractElement{Float64}}()
    loadset = Vector{NodeForce{Float64}}()

    #convert nodes
    for node in network.nodes
        pos = Vector(node.position)
        id = node.id

        #per-axis: FDM translation fixity carries over; rotations stay free
        #(truss rotations go inactive automatically)
        tn = Node(pos, vcat(collect(node.fixity), trues(3)))
        tn.id = id

        push!(nodeset, tn)
    end

    #convert elements
    for element in network.elements
        te = TrussElement(nodeset[[element.iStart, element.iEnd]]..., section)

        te.id = element.id

        push!(elset, te)
    end

    #convert loads
    for load in network.loads
        push!(loadset, NodeForce(nodeset[load.point.index], Vector(load.force)))
    end

    truss = Model(nodeset, elset, loadset)
    solve!(truss)

    return truss

end
