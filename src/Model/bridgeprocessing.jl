

function topologizematrix(elements::Vector{<:FrameElement}, indexer::Dict{Int64, Int64}, ordermatrix::SparseMatrixCSC{Float64, Int64})

    activeelements = Vector{Int64}()

    for (i, e) in enumerate(elements)
        rowstart = indexer[e.elementStart.elementID]
        rowend = indexer[e.elementEnd.elementID]

        push!(activeelements, rowstart, rowend)

        # ordermatrix[[rowstart, rowend], i + 1] .= [e.posStart, e.posEnd]
        ordermatrix[rowstart, i+1] = e.posStart
        ordermatrix[rowend, i+1] = e.posEnd
    end

    return sort!(unique(activeelements))
end

shatterReleaseDict = Dict(Element{FixedFixed} => (:fixedfixed, :fixedfixed),
    Element{FreeFixed} => (:freefixed, :fixedfixed),
    Element{FixedFree} => (:fixedfixed, :fixedfree),
    Element{FreeFixed} => (:freefixed, :fixedfree),
    Element{Joist} => (:freefixed, :fixedfree))

typeToRelease = Dict(BridgeElement{FixedFixed} => :fixedfixed,
    BridgeElement{FreeFixed} => :freefixed,
    BridgeElement{FixedFree} => :fixedfree,
    BridgeElement{FreeFree} => :freefree,
    BridgeElement{Joist} => :joist)

function shatter!(model::Model, 
        n::Int64,
        m::Int64,
        ordermatrix::SparseMatrixCSC{Float64, Int64}, 
        iactive::Vector{Int64}, 
        itrue::Vector{Int64},
        ibridge::Vector{Int64},
        indexerFwd::Dict{Int64, Int64},
        indexerRev::Dict{Int64, Int64}; tol = 1e-5)

    newnodes = Vector{Node}()
    newelements = Vector{Element}()

    nodematrix = spzeros(Int64, n, m)

    istart = model.nodes[end].nodeID

    i = 1

    #for each element that requires shattering
    for index in iactive

        #temp storage
        fracstore = Vector{Float64}()
        nodestore = Vector{Node}()
        elementstore = Vector{Element}()

        # element to shatter
        element = model.elements[indexerRev[index]]
        startpos = element.nodeStart.position
        # L = length(element)
        # vec = local_x(element)
        L = element.length
        vec = first(element.LCS)

        fractions = ordermatrix[index, 2:end-1]
        # istore = zeros(Int64, m)
        istore = Vector{Int64}()

        #iterate through fractures
        for (j, frac) in enumerate(fractions)
            frac == 0 && continue #no fracture occurs

            # if first non-zero value
            if isempty(fracstore)
                pos = startpos .+ vec * frac * L
                node = Node(pos, :free)
                node.nodeID = istart + i
                node.id = :bridge

                push!(nodestore, node)
                push!(fracstore, frac)
                push!(istore, i)

                nodematrix[index, j] = i
                # istore[j] = i

                i += 1
                continue
            end

            #for all other non-zero fracture values
            diffs = abs.(fracstore .- frac)
            val, ind = findmin(diffs)

            if val < tol
                nodematrix[index, j] = istore[ind] #node already exists
            else
                pos = startpos .+ vec * frac * L
                node = Node(pos, :free)
                node.nodeID = istart + i
                node.id = :bridge

                push!(nodestore, node)
                push!(fracstore, frac)
                push!(istore, i)

                nodematrix[index, j] = i
                # istore[j] = i

                i += 1
            end
        end

        # sort the new intermediate nodes along an element
        isorted = sortperm(fracstore)
        nodeset = [element.nodeStart; nodestore[isorted]; element.nodeEnd]


        releases = shatterReleaseDict[typeof(element)]

        # create intermediary elements
        for k = 1:length(nodeset) - 1
            if k == 1
                el = Element(nodeset[k], nodeset[k+1], element.section; release = releases[1])
            elseif k == length(nodeset) - 1
                el = Element(nodeset[k], nodeset[k+1], element.section; release = releases[2])
            else
                el = Element(nodeset[k], nodeset[k+1], element.section)
            end


            el.elementID = element.elementID
            el.Ψ = element.Ψ
            el.id = element.id

            push!(elementstore, el)
        end

        #collect all new nodes and new elements
        newnodes = [newnodes; nodestore]
        newelements = [newelements; elementstore]
    end


    #create new bridge elements
    for (k, be) in enumerate(model.elements[ibridge])
        rowstart = indexerFwd[be.elementStart.elementID]
        rowend = indexerFwd[be.elementEnd.elementID]

        nstart = newnodes[nodematrix[rowstart, k]]
        nend = newnodes[nodematrix[rowend, k]]

        el = Element(nstart, nend, be.section; release = typeToRelease[typeof(be)])
        el.elementID = be.elementID
        el.Ψ = be.Ψ
        el.id = be.id

        push!(newelements, el)
    end

    #delete shattered and bridge elements
    deleteat!(model.elements, sort!([itrue; ibridge]))

    model.elements = [model.elements; newelements]
    model.nodes = [model.nodes; newnodes]

    process_elements!(model)

end

function transfer!(load::LineLoad, id::Int64, inout::Bool, model::Model, elementids::Vector{Int64}, newloads::Vector{AbstractLoad})

    if !inout && typeof(load.element) <: Element
        return
    end

    itransfer = findall(elementids .== id)

    for i in itransfer
        newload = LineLoad(model.elements[i], load.value)
        newload.loadID = load.loadID
        newload.id = load.id
        push!(newloads, newload)
    end
end

function transfer!(load::PointLoad, id::Int64, inout::Bool, model::Model, elementids::Vector{Int64}, newloads::Vector{AbstractLoad})
    etype = typeof(load.element)
    typecheck = etype <: Element

    if !inout && typecheck
        return
    end


    if inout && typecheck #if load is applied to a shattered element
        position = load.element.length * load.position
        itransfer = findall(elementids .== id)

        endpositions = cumsum(getproperty.(model.elements[itransfer], :length))
        ielement = findfirst(endpositions .> position)

        newfrac = (position - endpositions[ielement-1]) / (endpositions[ielement] - endpositions[ielement-1])

        newload = PointLoad(model.elements[ielement], newfrac, load.value)
        newload.loadID = load.loadID
        newload.id = load.id

        push!(newloads, newload)
    else
        itransfer = findfirst(elementids .== id)

        newload = PointLoad(model.elements[itransfer], load.position, load.value)
        newload.loadID = load.loadID
        newload.id = load.id

        push!(newloads, newload)
    end
end


function convertloads!(model::Model, itrue::Vector{Int64})
    
    # collection of new loads
    newloads = Vector{AbstractLoad}()

    # indexes to remove
    rmid = Vector{Int64}()

    elementIDs = getproperty.(model.elements, :elementID)

    for (index, load) in enumerate(model.loads)

        #conversions only aplly to elemental loads
        typeof(load) <: NodeForce && continue

        #id of element load was applied to
        id = load.element.elementID

        #was this load applied to a now shattered element
        inout = in(id, itrue)

        #remove the existing load if applied to sh
        if inout || (typeof(load.element) <: BridgeElement)
            push!(rmid, index)
        end

        transfer!(load, id, inout, model, elementIDs, newloads)
    end

    #delete replaced loads
    deleteat!(model.loads, sort!(unique!(rmid)))

    model.loads = [model.loads; newloads]

end

function processBridge!(model::Model)

    #extract elements
    elements = model.elements

    #indices of regular and bridge elements
    iElement = findall(typeof.(elements) .<: Element)
    iBridge = findall(typeof.(elements) .<: BridgeElement)

    #get geometric information
    process_elements!(elements[iElement])

    #number of each type
    n = length(iElement)
    m = length(iBridge)

    #keep track of fractional shatter points for all element/bridge intersections
    ordermat = [spzeros(n, m+1) ones(n)]

    #convert index from all elements to index of Element type
    idvec = getproperty.(elements[iElement], :elementID)
    id2ind = Dict(idvec .=> collect(1:n))
    ind2id = Dict(collect(1:n) .=> idvec)

    #fill the order matrix and get indices of shattered elements
    iactive = topologizematrix(elements[iBridge], id2ind, ordermat)
    itrueactive = [ind2id[i] for i in iactive] #reverse index

    #shatter and make new nodes/elements
    shatter!(model, n, m, ordermat, iactive, itrueactive, iBridge, id2ind, ind2id; tol = model.tol)
    # println("Hello")

    #convert loads to new elements
    convertloads!(model, itrueactive)

    #generate DOF indices
    populate_DOF_indices!(model)

    #remake placeholders to reflect larger DOFs
    model.S = spzeros(Float64, 6model.nNodes, 6model.nNodes)
    model.P = zeros(Float64, 6model.nNodes)
    model.Pf = zeros(Float64, 6model.nNodes)
    model.u = zeros(Float64, 6model.nNodes)
    model.reactions = zeros(Float64, 6model.nNodes)

    # #generate global load vectors
    # populate_loads!(model)

    # #generate stiffnessmatrix
    # create_S!(model)

    # model.processed = true

end


# explicit version

"""
Test function for solving models with bridge elements
"""
function process_bridge!(model::Model)

    elements = model.elements

    iElement = findall(typeof.(elements) <: Element)
    iBridge = findall(typeof.(elements) .<: BridgeElement)

    process_elements!(Vector{Element}(elements[iElement]))

    n = length(iElement)
    m = length(iBridge)

    ordermat = [spzeros(n, m+1) ones(n)]

    idvec = getproperty.(elements[iElement], :elementID)

    id2ind = Dict(idvec .=> collect(1:n))
    ind2id = Dict(collect(1:n) .=> idvec)

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
    iactive = sort!(activeelements)
    itrueactive = [ind2id[i] for i in iactive]

    nodemat = spzeros(Int64, n, m)
    newnodes = Vector{Node}()
    newelements = Vector{Element}()

    nodeIDstart = model.nodes[end].nodeID

    i = 1
    tol = 1e-5
    for index in iactive
        fracstore = Vector{Float64}()
        nodestore = Vector{Node}()
        elementstore = Vector{Element}()

        #element to shatter
        element = elements[ind2id[index]]
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
    deleteat!(model.elements, sort!([itrueactive; iBridge])) # delete unshattered elements and all bridge elements
    model.elements = [model.elements; newelements]
    model.nodes = [model.nodes; newnodes]

    process_elements!(model)

    # process loads
    # convertloads!(model, itrueactive)
    newloads = Vector{Asap.AbstractLoad}()

    rmid = Vector{Int64}()

    elementids = getproperty.(model.elements, :elementID)
    for (index, load) in enumerate(model.loads)
        id = load.element.elementID
        in(id, itrueactive) && push!(rmid, index)

        if typeof(load) <: LineLoad
            # IF APPLIED TO A SHATTERED ELEMENT
            if in(id, itrueactive)
                itransfer = findall(elementids .== id)

                for i in itransfer
                    newload = LineLoad(model.elements[i], load.value)
                    newload.loadID = load.loadID
                    newload.id = load.id
                    push!(newloads, newload)
                end

            elseif typeof(load.element) <: BridgeElement #IF APPLIED TO A BRIDGE ELEMENT
                push!(rmid, index)

                itransfer = findfirst(elementids .== id)
                newload = LineLoad(model.elements[itransfer], load.value)
                newload.loadID = load.loadID
                newload.id = load.id
                push!(newloads, newload)
            end
        elseif typeof(load) <: PointLoad

            if typeof(load.element) <: Element

                if in(id, itrueactive)

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
    populate_DOF_indices!(model)

    model.S = spzeros(Float64, 6model.nNodes, 6model.nNodes)
    model.P = zeros(Float64, 6model.nNodes)
    model.Pf = zeros(Float64, 6model.nNodes)
    model.u = zeros(Float64, 6model.nNodes)
    model.reactions = zeros(Float64, 6model.nNodes)

    populate_loads!(model)

    Asap.create_S!(model)

    model.processed = true
end