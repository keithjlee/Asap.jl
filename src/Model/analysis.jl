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
function solve_explicit!(model::Model)

    elements = model.elements

    # find which elements are bridge elements
    iElement = findall(typeof.(elements) .== Element)
    iBridge = findall(typeof.(elements) == BridgeElement)

    # size of info arrays
    n = length(iElement)
    m = length(iBridge)

    # initialize arrays
    # C is the connectivity matrix, each column is a bridge element, -1: start element, 1: end element
    C = zeros(Int64, n, m)

    # fraction shatter point values for all bridge elements
    ordermat = [zeros(n, m+1) ones(n)]

    # get IDs of existing elements
    idvec = getproperty.(elements[iElement], :elementID)

    # convert element IDs to actual index in data matrices
    id2ind = Dict(idvec .=> collect(1:n))

    # populate data matrices for all bridge elements
    activeelements = Vector{Int64}()
    for (i, e) in enumerate(elements[BridgeElement])
        e1 = e.elementStart
        e2 = e.elementEnd

        # which elements (rows) does each bridge element connect to?
        rowstart = id2ind[e1.elementID]
        rowend = id2ind[e2.elementID]

        push!(iactive, rowstart, rowend)

        # populate connectivity matrix
        C[rowstart, i] = -1
        C[rowend, i] = 1

        # populate positional information
        ordermat[rowstart, i+1] = e.posStart
        ordermat[rowend, i+1] = e.posEnd
    end

    #relevant element indices
    unique!(activeelements)
    iactive = sortperm(activeelements)

    #generate all unique nodes and develop indices for shatter and bridge elements
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

            # if g is empty, frac value is unique
            if length(g) == 0 frac != 0
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
    deleteat!(model.elements, [idvec; iBridge]) # delete unshattered elements and all bridge elements
    model.elements = [model.elements; newelements]
    model.nodes = [model.nodes; newnodes]

    # process loads
    for load in model.loads
        if typeof(load) == LineLoad
            
        end

    end

end