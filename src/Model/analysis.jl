"""
    solve!(model::AbstractModel; reprocess = false)

Perform a structural analysis. `reprocess = true` re-generates the stiffness matrix and resets all saved solutions.
"""
function solve!(model::Model; reprocess = false)

    if !model.processed || reprocess
        process!(model)
    end

    # reduce scope of problem and solve
    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    U = model.S[idx, idx] \ F

    #compliance
    model.compliance = U' * F

    #full DOF displacement vector
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    # post process
    postprocess!(model)
end


function solve!(model::TrussModel; reprocess = false)

    if !model.processed || reprocess
        process!(model)
    end

    # reduce scope of problem and solve
    idx = model.freeDOFs
    U = model.S[idx, idx] \ model.P[idx]

    #compliance
    model.compliance = U' * model.P[idx]

    #full DOF displacement vector
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    # post process
    postprocess!(model)
end

"""
    solve(model::Model, L::Vector{Load})

Return the displacement vector to a new set of loads
"""
function solve(model::Model, L::Vector{<:Load})

    @assert model.processed "Model must be already processed"
    
    F = createF(model, L)
    
    idx = model.freeDOFs

    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
    solve!(model::Model, L::Vector{Load})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::Model, L::Vector{<:Load})

    @assert model.processed "Model must be already processed"

    model.loads = L

    # clear existing load associations
    for (node, element) in zip(model.nodes, model.elements)
        empty!(element.loadIDs)
    end

    # assign new load associations
    for (i, load) in enumerate(model.loads)
        load.loadID = i
        assign!(load)
    end

    # create P, Pf matrix
    populateLoads!(model)
    
    # reduce scope of problem and solve
    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    U = model.S[idx, idx] \ F

    #compliance
    model.compliance = U' * F

    #full DOF displacement vector
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    # post process
    postprocess!(model)
end

"""
    solve(model::TrussModel, L::Vector{NodeForce})

Return the displacement vector to a new set of loads
"""
function solve(model::TrussModel, L::Vector{NodeForce})

    !model.processed && process!(model)
    
    F = createF(model, L)
    
    idx = model.freeDOFs

    U = model.S[idx, idx] \ F[idx]

    u = zeros(model.nDOFs)
    u[idx] = U

    return u
end

"""
    solve!(model::TrussModel, L::Vector{NodeForce})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::TrussModel, L::Vector{NodeForce})

    @assert model.processed "Model must be already processed"

    model.loads = L

    # clear existing load associations
    for element in model.elements
        empty!(element.loadIDs)
    end

    # assign new load associations
    for (i, load) in enumerate(model.loads)
        load.loadID = i
        assign!(load)
    end

    # create P, Pf matrix
    populateLoads!(model)
    
    # reduce scope of problem and solve
    idx = model.freeDOFs
    U = model.S[idx, idx] \ model.P[idx]

    #compliance
    model.compliance = U' * model.P[idx]

    #full DOF displacement vector
    model.u = zeros(model.nDOFs)
    model.u[idx] = U

    # post process
    postprocess!(model)
end

"""
Test function for solving models with bridge elements
"""
function process_bridge!(model::Model)

    elements = model.elements

    iElement = findall(typeof.(elements) .== Element)
    iBridge = findall(typeof.(elements) .== BridgeElement)

    processElements!(Vector{Element}(elements[iElement]))

    n = length(iElement)
    m = length(iBridge)

    ordermat = [zeros(n, m+1) ones(n)]

    idvec = getproperty.(elements[iElement], :elementID)

    id2ind = Dict(idvec .=> collect(1:n))

    activeelements = Vector{Int64}()
    for (i, e) in enumerate(elements[iBridge])
        e1 = e.elementStart
        e2 = e.elementEnd

        # which elements (rows) does each bridge element connect to?
        rowstart = id2ind[e1.elementID]
        rowend = id2ind[e2.elementID]

        push!(activeelements, rowstart, rowend)

        # populate positional information
        ordermat[rowstart, i+1] = e.posStart
        ordermat[rowend, i+1] = e.posEnd
    end

    #relevant element indices
    unique!(activeelements)
    iactive = sortperm(activeelements)

    nodemat = zeros(Int64, n, m)
    newnodes = Vector{Node}()
    newelements = Vector{Element}()

    nodeIDstart = model.nodes[end].nodeID

    i = 1
    tol = 1e-2
    for index in iactive
        fracstore = Vector{Float64}()
        nodestore = Vector{Node}()
        elementstore = Vector{Element}()

        #element to shatter
        element = elements[index]
        startpos = element.nodeStart.position
        L = element.length
        vec = first(element.LCS)

        fractions = ordermat[index, 2:end-1]
        istore = zeros(Int64, length(fractions))

        for (j, frac) in enumerate(fractions)

            frac == 0 && continue

            # if fracstore is empty, frac value is unique
            if length(fracstore) == 0
                pos = startpos .+ vec * frac * L
                node = Node(pos, :free)
                node.nodeID = nodeIDstart + i
                node.id = :bridge

                push!(nodestore, node)
                push!(fracstore, frac)

                nodemat[index, j] = i
                istore[j] = i
                
                i += 1
                continue
            end

            # find the difference between this fracture value and previous values
            diffs = abs.(fracstore .- frac)

            # find minimum difference
            val, ind = findmin(diffs)

            # if min distance is less than tolerance, it already exists as a node
            if val < tol
                nodemat[index, j] = istore[ind]
            else # create a new node
                pos = startpos .+ vec * frac * L
                node = Node(pos, :free)
                node.nodeID = nodeIDstart + i
                node.id = :bridge

                push!(nodestore, node)
                push!(fracstore, frac)

                nodemat[index, j] = i
                istore[j] = i
                
                i += 1
            end

        end

        # make shattered elements
        isorted = sortperm(fracstore)
        nodeset = [element.nodeStart; nodestore[isorted]; element.nodeEnd]

        for i = 1:length(nodeset) - 1
            el = Element(nodeset[i], nodeset[i+1], element.section)
            el.elementID = element.elementID
            el.Ψ = element.Ψ
            el.id = element.id

            push!(elementstore, el)
        end

        # collect all new elements and nodes
        newnodes = [newnodes; nodestore]
        newelements = [newelements; elementstore]

    end

    # create bridge elements
    for (i, be) in enumerate(elements[iBridge])
        rowstart = id2ind[be.elementStart.elementID]
        rowend = id2ind[be.elementEnd.elementID]

        nstart = newnodes[nodemat[rowstart, i]]
        nend = newnodes[nodemat[rowend, i]]

        el = Element(nstart, nend, be.section)
        el.elementID = be.elementID
        el.Ψ = be.Ψ
        el.id = be.id
        el.release = be.release

        push!(newelements, el)
    end

    # modify input model with new elements and nodes
    deleteat!(model.elements, sort!([iactive; iBridge])) # delete unshattered elements and all bridge elements
    model.elements = [model.elements; newelements]
    model.nodes = [model.nodes; newnodes]

    processElements!(model)

    # process loads
    newloads = Vector{Asap.Load}()

    rmid = Vector{Int64}()

    elementids = getproperty.(model.elements, :elementID)
    for (index, load) in enumerate(model.loads)
        id = load.element.elementID
        if typeof(load) == LineLoad
            # IF APPLIED TO A SHATTERED ELEMENT
            if in(id, iactive)
                push!(rmid, index)
                itransfer = findall(elementids .== id)

                for i in itransfer
                    newload = LineLoad(model.elements[i], load.value)
                    newload.loadID = load.loadID
                    newload.id = load.id
                    push!(newloads, newload)
                end

            elseif typeof(load.element) == BridgeElement #IF APPLIED TO A BRIDGE ELEMENT
                push!(rmid, index)

                itransfer = findfirst(elementids .== id)
                newload = LineLoad(model.elements[itransfer], load.value)
                newload.loadID = load.loadID
                newload.id = load.id
                push!(newloads, newload)
            end
        elseif typeof(load) == PointLoad

            if typeof(load.element) == Element

                if in(id, iactive)
                    push!(rmid, index)

                    position = load.element.length * load.position

                    itransfer = findall(elementids .== id)

                    endpositions = cumsum(getproperty.(model.elements[itransfer], :length))
                    ielement = findfirst(endpositions .> position)
                    newfrac = (position - endpositions[ielement-1]) / (endpositions[ielement] - endpositions[ielement - 1])

                    newload = PointLoad(model.elements[itransfer[ielement]], newfrac, load.value)
                    newload.loadID = load.loadID
                    newload.id = load.id
                    push!(newloads, newload)
                end

            else
                push!(rmid, index)
                itransfer = findfirst(elementids .== id)
                newload = PointLoad(model.elements[itransfer], load.position, load.value)
                newload.loadID = load.loadID
                newload.id = load.id
                push!(newloads, newload)
            end

        end

    end

    deleteat!(model.loads, sort!(unique!(rmid)))
    model.loads = [model.loads; newloads]


    # populating DOFs
    model.DOFs = vcat([node.dof for node in model.nodes]...)
    model.nDOFs = length(model.DOFs)
    model.nNodes = length(model.nodes)
    model.nElements = length(model.elements)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)

    n_dof = 6
    dofset = collect(0:n_dof-  1)

    # assign an id to node, extract global DOF index
    @inbounds for (i, node) in enumerate(model.nodes)
        node.globalID = i * n_dof - (n_dof - 1) .+ dofset
    end

    #assign an id to load, store load id into relevant node/element
    @inbounds for load in model.loads
        Asap.assign!(load)
    end

    # assign an id to element, get associated node IDs, extract global DOF
    @inbounds for element in model.elements
        element.nodeIDs = [element.nodeStart.nodeID, element.nodeEnd.nodeID]

        idStart = element.nodeStart.globalID
        idEnd = element.nodeEnd.globalID

        element.globalID = [idStart; idEnd]
    end

    populateLoads!(model)

    Asap.globalS!(model)

    model.processed = true
end