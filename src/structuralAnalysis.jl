# global axes
const globalX = [1., 0., 0.]
const globalY = [0., 1., 0.]
const globalZ = [0., 0., 1.]

# expands elemental start/end nodes to global DOF indices
function dofExpander!(element::Element, nodes::Vector{Node})
    element.dofIndex = vcat(nodes[element.nodeIndex[1]].globalIndex, nodes[element.nodeIndex[2]].globalIndex)
end

# dofExpander for a vector of elements
function dofExpander!(elements::Vector{Element}, nodes::Vector{Node})
    for element in elements
        element.dofIndex = vcat(nodes[element.nodeIndex[1]].globalIndex, nodes[element.nodeIndex[2]].globalIndex)
    end
end

# global dof index of each node
function nodeGlobalIndex!(nodes::Vector{Node})
    n = length(nodes)
    nNodalDOFS = length(nodes[1].DOFS)
    for i = 1:n
        nodes[i].globalIndex = i * nNodalDOFS - (nNodalDOFS-1) .+ collect(0:nNodalDOFS-1)
    end
end



```
Local coordinate system of element
```
function lcs(element::Element, Ψ; tol = 0.001)

    # local x vector
    xvec = normalize(element.posEnd .- element.posStart)
    
    if norm(cross(xvec, globalY)) < tol
        CYx = xvec[2] #cosine to global Y axis
        xvec = CYx * globalY
        yvec = -CYx * globalX * cos(Ψ) + sin(Ψ) * globalZ
        zvec = CYx * globalX * sin(Ψ) + cos(Ψ) * globalZ 
    else
        zbar = normalize(cross(xvec, [0, 1, 0]))
        ybar = normalize(cross(zbar, xvec))

        yvec = cos(Ψ) * ybar + sin(Ψ) * zbar
        zvec = -sin(Ψ) * ybar + cos(Ψ) * zbar
    end

    return Vec3.([xvec, yvec, zvec])
end


#scales displacement s.t. maximum displacement is 30% of the mean element length
function autoScaleFactor(structure::Structure)
    lengths = [e.length for e in structure.elements]
    if structure.dims == 2
        displacements = vcat([structure.U[n.globalIndex[1:2]] for n in structure.nodes]...)
    else
        displacements = vcat([structure.U[n.globalIndex[1:3]] for n in structure.nodes]...)
    end
    factor = 0.30 / maximum(abs.(displacements)) * mean(lengths)
    return round(factor)
end

function postProcess!(structure::Structure; scaleFactor = :auto)

    ### create displaced nodes and elements (scaled)
    if scaleFactor == :auto
        scaleFactor = autoScaleFactor(structure)
    end

    positions = [node.position for node in structure.nodes]

    if structure.dims == 2
        displacements = vcat([structure.U[n.globalIndex[1:2]] for n in structure.nodes])
    else
        displacements = vcat([structure.U[n.globalIndex[1:3]] for n in structure.nodes])
    end

    newpositions = [positions[i] .+ scaleFactor .* displacements[i] for i = 1:length(structure.nodes)]

    structure.displacedNodes = Node.(newpositions)
    structure.displacedElements = [Element(structure.displacedNodes, element.nodeIndex, element.type) for element in structure.elements]

    ### Internal forces in local CS
    for element in structure.elements
        element.localForce = element.R * element.k * structure.U[element.dofIndex]
        element.globalForce = element.k * structure.U[element.dofIndex]
        element.axialForce = axialLoad(element)
    end

    ### Nodal displacements and reactions
    for node in structure.nodes
        node.disp = structure.U[node.globalIndex]
        if any(node.DOFS .== false)
            node.reaction = structure.reactions[node.globalIndex]
        end
    end
end

function reactions!(structure::Structure)
    if isdefined(structure, :U) == false
        error("Perform analysis before solving reactions.")
    end

    # index of boundary condition DOFs
    fixedDOFS = findall(structure.DOFS .== false)

    # solve for loads at reactions
    reactions = structure.K[fixedDOFS, :] * structure.U

    # initialize
    structure.reactions = zeros(structure.nDOFS)
    
    # update
    structure.reactions[fixedDOFS] = reactions
end


function addNodeElements!(elements::Vector{Element}, nodes::Vector{Node})
    for i = 1:length(elements)
        j, k = elements[i].nodeIndex
        push!(nodes[j].elements, (i, -1))
        push!(nodes[k].elements, (i, 1))
    end

    if isdefined(nodes[1], :elements)
        for i = 1:length(nodes)
            nodes[i].elements = unique(nodes[i].elements)
        end
    end
end

function addNodeLoads!(loads::Vector{Load}, nodes::Vector{Node})
   for i = 1:length(loads)
    nodes[loads[i].index].load = loads[i].load 
   end
end

function analyze!(structure::Structure; forceK = false, SF = :auto)
    if !isdefined(structure, :F)
        error("No loads defined.")
    end

    # associated elements to each node
    addNodeElements!(structure.elements, structure.nodes)
    # associated load vectors to individual nodes
    addNodeLoads!(structure.loads, structure.nodes)

    # prevents rebuilding of stiffness matrix if structure stays the same
    # setting forceK = true rebuilds the stiffness matrix
    if isdefined(structure, :K) && !forceK
        #displacement
        U = structure.K[structure.freeDOFS, structure.freeDOFS] \ structure.F[structure.freeDOFS]

        #compliance
        structure.compliance = U' * structure.F[structure.freeDOFS]
        return
    end
    #create nodal DOF indices
    nodeGlobalIndex!(structure.nodes)

    #create elemental DOF indices
    dofExpander!(structure.elements, structure.nodes)

    #create elemental stiffness matrices
    k_elemental!.(structure.elements, structure.dims)

    #create global stiffness matrix
    structure.K = K(structure.elements, structure.nDOFS)

    #displacement
    U = structure.K[structure.freeDOFS, structure.freeDOFS] \ structure.F[structure.freeDOFS]

    #compliance
    structure.compliance = U' * structure.F[structure.freeDOFS]

    #long form displacement
    structure.U = zeros(structure.nDOFS)
    structure.U[structure.freeDOFS] = U

    #reactions
    reactions!(structure)

    #post processing
    postProcess!(structure; scaleFactor = SF)
end


function axialLoad(element::Element)
    if isdefined(element, :globalForce) == false
        error("No internal loads defined; run analysis and/or postprocessing")
    end

    dims = length(element.posStart)

    if element.type == :truss
        if dims == 2
            if dot(element.posEnd - element.posStart, element.globalForce[3:4]) < 0
                return - norm(element.globalForce[3:4])
            else
                return norm(element.globalForce[3:4])
            end
        else
            if dot(element.posEnd - element.posStart, element.globalForce[4:6]) < 0
                return - norm(element.globalForce[4:6])
            else
                return norm(element.globalForce[4:6])
            end
        end
    else
        if dims == 2
            if dot(element.posEnd - element.posStart, element.globalForce[4:5]) < 0
                return - norm(element.globalForce[4:5])
            else
                return norm(element.globalForce[4:5])
            end
        else
            if dot(element.posEnd - element.posStart, element.globalForce[7:9]) < 0
                return - norm(element.globalForce[7:9])
            else
                return norm(element.globalForce[7:9])
            end
        end
    end
end

function trussSizer(P, k, L, E, ρ, σ, minA)
    
    Astress = P/σ

    if sign(P) < 0
        D = (64 * abs(P) * (k * L)^2 / pi^3 / E / (1-ρ^4))^(1/4)
        Abuckling = pi/4 * D^2 *(1-ρ^2)
        if !isnothing(minA)
            A = max(Abuckling, Astress, minA)
        else
            A = max(Abuckling, Astress)
        end
    else
        if !isnothing(minA)
            A = max(Astress, minA)
        else
            A = Astress
        end
    end

    return A
end

# pseudo size is an approximated optimization function for sizing member areas. It performs a linear elastic analysis using uniform section properties, then resizes members to meet stress and buckling limitiations
function pseudoSize!(structure::Structure,
        maxStress;
        minArea = nothing,
        dDratio = 0.9,
        k = 1,
        maxIters = 1)

    # perform analysis

    analyze!(structure)

    iter = 0
    while iter < maxIters
        for element in structure.elements
            element.A = trussSizer(element.axialForce, k, element.length, element.E, dDratio, maxStress, minArea)
        end
        analyze!(structure; forceK = true)
        iter += 1
    end
end

function pseudoSize(structure::Structure,
    maxStress;
    minArea = nothing,
    dDratio = 0.9,
    k = 1,
    maxIters = 1)

    # perform analysis
    newstructure = deepcopy(structure)
    analyze!(newstructure)

    iter = 0
    while iter < maxIters
        if storeAreas
            areas = Vector{Float64}()
        end
        for element in newstructure.elements
            element.A = trussSizer(element.axialForce, k, element.length, element.E, dDratio, maxStress, minArea)
            if storeAreas
                push!(areas, element.A)
            end
        end
        analyze!(newstructure; forceK = true)

        iter += 1
    end

    return newstructure
end

function structureMass(structure::Structure, ρ)
    volume = 0.0
    for element in structure.elements
        volume += element.A * element.length
    end
    return volume * ρ
end