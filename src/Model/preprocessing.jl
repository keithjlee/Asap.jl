"""
populate global DOF indices of nodes and elements
"""
function populateDOF!(model::Model)

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
        assign!(load)
    end

    # assign an id to element, get associated node IDs, extract global DOF
    @inbounds for element in model.elements
        element.nodeIDs = [element.nodeStart.nodeID, element.nodeEnd.nodeID]

        idStart = element.nodeStart.globalID
        idEnd = element.nodeEnd.globalID

        element.globalID = [idStart; idEnd]
    end

end

"""
populate global DOF indices of nodes and elements
"""
function populateDOF!(model::TrussModel)

    model.DOFs = vcat([node.dof for node in model.nodes]...)
    model.nDOFs = length(model.DOFs)
    model.nNodes = length(model.nodes)
    model.nElements = length(model.elements)
    model.freeDOFs = findall(model.DOFs)
    model.fixedDOFs = findall(.!model.DOFs)

    n_dof = 3
    @inbounds for (i, node) in enumerate(model.nodes)
        node.globalID = i * n_dof - (n_dof - 1) .+ collect(0:n_dof - 1)
    end

    @inbounds for load in model.loads
        assign!(load)
    end

    @inbounds for element in model.elements

        element.nodeIDs = [element.nodeStart.nodeID, element.nodeEnd.nodeID]

        idStart = element.nodeStart.globalID
        idEnd = element.nodeEnd.globalID

        element.globalID = [idStart; idEnd]
    end

end

"""
Process elements: get transformation matrix and global elemental stiffness matrix
"""
function processElements!(model::Model)
    for element in model.elements
        element.Q = zeros(12) # reset Qf
        element.R = R(element)
        element.LCS = lcs(element, element.Ψ)
        element.length = length(element)
        makeK!(element)
    end
end

function processElements!(elements::Vector{<:FrameElement})
    for element in elements
        element.Q = zeros(12) # reset Qf
        element.R = R(element)
        element.LCS = lcs(element, element.Ψ)
        element.length = length(element)
        makeK!(element)
    end
end

"""
Process elements: get transformation matrix and global elemental stiffness matrix
"""
function processElements!(model::TrussModel)
    for element in model.elements
        element.R = R(element)
        element.LCS = lcs(element, element.Ψ)
        element.length = length(element)
        makeK!(element)
    end
end

"""
Populate a node force load
"""
function populateLoad!(model::Model, load::NodeForce)
    idx = load.node.globalID[1:3]
    model.P[idx] += load.value
end

"""
Populate a node force load
"""
function populateLoad!(model::TrussModel, load::NodeForce)
    idx = load.node.globalID
    model.P[idx] += load.value
end


function populateLoad!(P::Vector{Float64}, load::NodeForce)
    idx = load.node.globalID[1:3]
    P[idx] += load.value
end

"""
Populate a node moment load
"""
function populateLoad!(model::Model, load::NodeMoment)
    idx = load.node.globalID[4:6]
    model.P[idx] += load.value
end

function populateLoad!(P::Vector{Float64}, load::NodeMoment)
    idx = load.node.globalID[4:6]
    P[idx] += load.value
end

"""
Populate element loads
"""
function populateLoad!(model::Model, load::ElementLoad)
    idx = load.element.globalID
    load.element.Q += load.element.R' * Q(load)
    model.Pf[idx] += load.element.Q
end

function populateLoad!(P::Vector{Float64}, load::ElementLoad)
    idx = load.element.globalID
    P[idx] += load.element.R' * Q(load)
end

"""
Create load vectors P and Pf
"""
function populateLoads!(model::Model)
    #initialize
    model.P = zeros(model.nDOFs)
    model.Pf = zeros(model.nDOFs)

    #create load vectors
    for load in model.loads
        populateLoad!(model, load)
    end
end

"""
Create load vectors P and Pf
"""
function populateLoads!(model::TrussModel)
    #initialize
    model.P = zeros(model.nDOFs)

    #create load vectors
    for load in model.loads
        populateLoad!(model, load)
    end
end

"""
create load vector F = P-Pf
"""
function createF(model::Model, loads::Vector{Load})

    P = zeros(model.nDOFs)
    Pf = zeros(model.nDOFs)

    for load in loads
        if typeof(load) <: ElementLoad
            populateLoad!(Pf, load)
        else
            populateLoad!(P, load)
        end
    end

    return P - Pf
end

"""
create load vector F = P-Pf
"""
function createF(model::TrussModel, loads::Vector{NodeForce})

    P = zeros(model.nDOFs)

    for load in loads
        populateLoad!(P, load)
    end

    return P
end

"""
Create global stiffness matrix
"""
function globalS!(model::Model)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()

    for element in model.elements

        idx = element.globalID

        @inbounds for i = 1:12
            @inbounds for j = 1:12
                k = element.K[i,j]
                push!(I, idx[i])
                push!(J, idx[j])
                push!(V, k)
                
            end
        end

    end

    model.S = sparse(I, J, V)
end

"""
Create global stiffness matrix
"""
function globalS!(model::TrussModel)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()

    for element in model.elements

        idx = element.globalID

        @inbounds for i = 1:6
            @inbounds for j = 1:6
                k = element.K[i,j]
                push!(I, idx[i])
                push!(J, idx[j])
                push!(V, k)
                
            end
        end

    end

    model.S = sparse(I, J, V)
end

"""
process a network
"""
function process!(model::AbstractModel)

    #elements 
    processElements!(model)

    #global DOF 
    populateDOF!(model)

    #loads
    populateLoads!(model)

    #stiffness matrix
    globalS!(model)

    #processing finished
    model.processed = true
end



function topologizematrix(elements::Vector{BridgeElement}, indexer::Dict{Int64, Int64}, ordermatrix::SparseMatrixCSC{Float64, Int64})

    activeelements = Vector{Int64}()

    @inbounds for (i, e) in enumerate(elements)
        rowstart = indexer[e.elementStart.elementID]
        rowend = indexer[e.elementEnd.elementID]

        push!(activeelements, rowstart, rowend)

        ordermatrix[[rowstart, rowend], i + 1] .= [e.posStart, e.posEnd]
    end

    return sort!(unique(activeelements))
end

shatterReleaseDict = Dict(:fixedfixed => (:fixedfixed, :fixedfixed),
    :freefixed => (:freefixed, :fixedfixed),
    :fixedfree => (:fixedfixed, :fixedfree),
    :freefree => (:freefixed, :fixedfree))

function shatter!(model::Model, 
        n::Int64,
        m::Int64,
        ordermatrix::SparseMatrixCSC{Float64, Int64}, 
        iactive::Vector{Int64}, 
        itrue::Vector{Int64},
        ibridge::Vector{Int64},
        indexerFwd::Dict{Int64, Int64},
        indexerRev::Dict{Int64, Int64}; tol = 1e-2)

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
        L = element.length
        vec = first(element.LCS)

        fractions = ordermatrix[index, 2:end-1]
        istore = zeros(Int64, m)

        #iterate through fractures
        for (j, frac) in enumerate(fractions)
            frac == 0 && continue #no fracture occurs

            # if first non-zero value
            if length(fracstore) == 0
                pos = startpos .+ vec * frac * L
                node = Node(pos, :free)
                node.nodeID = istart + 1
                node.id = :bridge

                push!(nodestore, node)
                push!(fracstore, frac)

                nodematrix[index, j] = i
                istore[j] = i

                i += 1
            end

            #for all other non-zero fracture values
            diffs = abs.(fracstore .- frac)
            val, ind = findmin(diffs)

            if val < tol
                nodematrix[index, j] = istore[ind] #node already exisdts
            else
                pos = startpos .+ vec * frac * L
                node = Node(pos, :free)
                node.nodeID = istart + 1
                node.id = :bridge

                push!(nodestore, node)
                push!(fracstore, frac)

                nodematrix[index, j] = i
                istore[j] = i

                i += 1
            end
        end

        # sort the new intermediate nodes along an element
        isorted = sortperm(fracstore)
        nodeset = [element.nodeStart; nodestore[isorted]; element.nodeEnd]


        releases = shatterReleaseDict[element.release]

        # create intermediary elements
        for i = 1:length(nodeset) - 1
            el = Element(nodeset[i], nodeset[i+1], element.section)
            el.elementID = element.elementID
            el.Ψ = element.Ψ
            el.id = element.id

            if i == 1
                el.release = releases[1]
            elseif i == length(nodeset) - 1
                el.release = releases[2]
            end

            push!(elementstore, el)
        end

        #collect all new nodes and new elements
        newnodes = [newnodes; nodestore]
        newelements = [newelements; elementstore]
    end


    #create new bridge elements
    for (i, be) in enumerate(model.elements[ibridge])
        rowstart = indexerFwd[be.elementStart.elementID]
        rowend = indexerFwd[be.elementEnd.elementID]

        nstart = newnodes[nodematrix[rowstart, i]]
        nend = newnodes[nodematrix[rowend, i]]

        el = Element(nstart, nend, be.section)
        el.elementID = be.elementID
        el.Ψ = be.Ψ
        el.id = be.id
        el.release = be.release

        push!(newelements, el)
    end

    #delete shattered and bridge elements
    deleteat!(model.elements, sort!([itrue; ibridge]))

    model.elements = [model.elements; newelements]
    model.nodes = [model.nodes; newnodes]

    processElements!(model)

end

function transfer!(load::LineLoad, id::Int64, inout::Bool, model::Model, elementids::Vector{Int64}, newloads::Vector{Load})

    if !inout && typeof(load.element) == Element
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

function transfer!(load::PointLoad, id::Int64, inout::Bool, model::Model, elementids::Vector{Int64}, newloads::Vector{Load})
    etype = typeof(load.element)
    typecheck = etype == Element

    if !inout && typecheck
        return
    end

    if typecheck
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


function convertloads!(model::Model, iactive::Vector{Int64})
    
    # collection of new loads
    newloads = Vector{Load}()

    # indexes to remove
    rmid = Vector{Int64}()

    elementIDs = getproperty.(model.elements, :elementID)

    for (index, load) in enumerate(model.loads)

        #conversions only aplly to elemental loads
        typeof(load) <: NodeForce && continue

        #id of element load was applied to
        id = load.element.elementID

        #was this load applied to a now shattered element
        inout = in(id, iactive)

        #remove existing load if it is applied on a shattered element or a bridge element
        inout || typeof(load.element) == BridgeElement && push!(rmid, index)

        transfer!(load, id, inout, model, elementIDs, newloads)
    end

    #delete replaced loads
    deleteat!(model.loads, sort!(unique!(rmid)))

    model.loads = [model.loads; newloads]

end

function processbridge(model::Model)

    #extract elements
    elements = model.elements

    #indices of regular and bridge elements
    iElement = findall(typeof.(elements) .== Element)
    iBridge = findall(typeof.(elements) .== BridgeElement)

    #get geometric information
    processElements!(Vector{Element}(elements[iElement]))

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
    iactive = topologizematrix(Vector{BridgeElement}(elements[iBridge]), id2ind, ordermat)
    itrueactive = [ind2id[i] for i in iactive] #reverse index

    #shatter and make new nodes/elements
    shatter!(model, n, m, ordermat, iactive, itrueactive, iBridge, id2ind, ind2id)
    
    #convert loads to new elements
    convertloads!(model, itrueactive)

    #generate DOF indices
    populateDOF!(model)

    #generate global load vectors
    populateLoads!(model)

    #generate stiffnessmatrix
    globalS!(model)

    model.processed = true

end