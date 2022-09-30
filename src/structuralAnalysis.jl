# global axes
const globalX = [1., 0., 0.]
const globalY = [0., 1., 0.]
const globalZ = [0., 0., 1.]

# expands elemental start/end nodes to global DOF indices
function dofExpander(element::Element, nodes::Vector{Node})
    element.dofIndex = vcat(nodes[element.nodeIndex[1]].globalIndex, nodes[element.nodeIndex[2]].globalIndex)
end

# dofExpander for a vector of elements
function dofExpander(elements::Vector{Element}, nodes::Vector{Node})
    for element in elements
        element.dofIndex = vcat(nodes[element.nodeIndex[1]].globalIndex, nodes[element.nodeIndex[2]].globalIndex)
    end
end

# global dof index of each node
function nodeGlobalIndex(nodes::Vector{Node})
    n = length(nodes)
    nNodalDOFS = length(nodes[1].DOFS)
    for i = 1:n
        nodes[i].globalIndex = i * nNodalDOFS - (nNodalDOFS-1) .+ collect(0:nNodalDOFS-1)
    end
end

```
SHAPE FUNCTIONS gives displacement orthogonal to local x axis based on end moments and shear forces
u1 = shear displacement at beginning node
u2 = moment displacement at beginning node
u3 = shear displacement at end node
u4 = moment displacement at end node

These displacements should be w/r/t local coordinate system
```

# 
function N1(x, element::Element)
    return 1 - 3 * (x / element.length)^2 + 2 * (x / element.length)^3
end

function N2(x, element::Element)
    return x * (1 - x / element.length)^2
end

function N3(x, element::Element)
    return 3 * (x / element.length)^2 - 2 * (x / element.length)^3
end

function N4(x, element::Element)
    return x^2 / element.length * (x / element.length - 1)
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

# creates the elemental stiffness matrix in global coordinate system
function k_elemental(element::Element, dims::Int; return_k = false, tol = 1e-3)
    type = element.type
    if type == :truss
        if dims == 3
            Cx, Cy, Cz = (element.posEnd .- element.posStart) ./ element.length

            # rotation matrix
            R = [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]
            element.R = R
            
            # stiffness matrix
            K = element.E * element.A / element.length .* [1 -1; -1 1]

            element.k = Symmetric(R' * K * R)
        elseif dims == 2
            Cx, Cy = (element.posEnd .- element.posStart) ./ element.length

            # rotation matrix
            R = [Cx Cy 0 0; 0 0 Cx Cy]
            element.R = R

            # stiffness matrix
            K = element.E * element.A / element.length .* [1 -1; -1 1]

            element.k = Symmetric(R' * K * R)
        else
            error("dims must be 2 or 3")
        end
    elseif type == :frame
        if dims == 3

            ##
            # deprecated
            ##
            #rotation matrix
            # CXx, CYx, CZx = (element.posEnd .- element.posStart) ./ element.length
            # D = sqrt(CXx^2 + CYx^2)
            # CXy = -CYx / D
            # CYy = CXx / D
            # CZy = 0
            # CXz = -CXx * CZx / D
            # CYz = -CYx * CZx / D
            # CZz = D
            
            # if abs(element.posStart[1]-element.posEnd[1]) < tol && abs(element.posStart[2]-element.posEnd[2]) < tol
            #     if element.posEnd[3] > element.posStart[3]
            #         Λ = [0 0 1; 0 1 0; -1 0 0]
            #     else
            #         Λ = [0 0 -1; 0 1 0; 1 0 0]
            #     end
            # else
            #     Λ = [CXx CYx CZx; CXy CYy CZy; CXz CYz CZz]
            # end

            ##
            # New method provides precise pitch angle control w/r/t local x axis
            xvec = normalize(element.posEnd .- element.posStart) # local x vector
            CXx, CYx, CZx = xvec # local x cosines

            
            if norm(cross(xvec, globalY)) < tol #special case for horizontal members aligned with global Y
                Λ = [0. CYx 0.;
                    -CYx*cos(element.Ψ) 0 sin(element.Ψ);
                    CYx*sin(element.Ψ) 0 cos(element.Ψ)]
            else # all other
                b1 = (-CXx * CYx * cos(element.Ψ) - CZx * sin(element.Ψ)) / sqrt(CXx^2 + CZx^2)
                b2 = sqrt(CXx^2 + CZx^2) * cos(element.Ψ)
                b3 = (-CYx * CZx * cos(element.Ψ) + CXx * sin(element.Ψ)) / sqrt(CXx^2 + CZx^2)

                c1 = (CXx * CYx * sin(element.Ψ) - CZx * cos(element.Ψ)) / sqrt(CXx^2 + CZx^2)
                c2 = -sqrt(CXx^2 + CZx^2) * sin(element.Ψ)
                c3 = (CYx * CZx * sin(element.Ψ) + CXx * cos(element.Ψ)) / sqrt(CXx^2 + CZx^2)

                Λ = [CXx CYx CZx; 
                    b1 b2 b3; 
                    c1 c2 c3]
            end
            
            R = [Λ zeros(3,9); zeros(3,3) Λ zeros(3,6); zeros(3,6) Λ zeros(3,3); zeros(3,9) Λ]
            element.R = R

            # Faulty method from Ferreira "MATLAB Codes for Finite Element Analysis"
            #stiffness matrix
            #nodal stiffnesses
            # k1 = element.E * element.A / element.length
            # k2 = 12 * element.E * element.Iz / element.length^3
            # k3 = 6 * element.E * element.Iz / element.length^2
            # k4 = 4 * element.E * element.Iz / element.length
            # k5 = 2 * element.E * element.Iz / element.length
            # k6 = 12 * element.E * element.Iy / element.length^3
            # k7 = 6 * element.E * element.Iy / element.length^2
            # k8 = 4 * element.E * element.Iy / element.length
            # k9 = 2 * element.E * element.Iy / element.length
            # k10 = element.G * element.J / element.length

            # #stiffness blocks
            # a = [k1 0 0; 0 k2 0; 0 0 k6]
            # b = [0 0 0; 0 0 k3; 0 -k7 0]
            # c = [k10 0 0; 0 k8 0; 0 0 k4]
            # d = [-k10 0 0; 0 k9 0; 0 0 k5]

            # element.k = Symmetric(R' * [a b -a b; b' c b d; -a' b' a -b; b' d' -b' c] * R)

            E = element.E
            length = element.length
            A = element.A
            Iz = element.Iz
            Iy = element.Iy
            G = element.G
            J = element.J
            k = E / length^3 * [
                A*length^2 0 0 0 0 0 -A*length^2 0 0 0 0 0;
                0 12Iz 0 0 0 6length*Iz 0 -12Iz 0 0 0 6length*Iz;
                0 0 12Iy 0 -6length*Iy 0 0 0 -12Iy 0 -6length*Iy 0;
                0 0 0 G*J*length^2/E 0 0 0 0 0 -G*J*length^2/E 0 0;
                0 0 -6length*Iy 0 4length^2*Iy 0 0 0 6length*Iy 0 2length^2*Iy 0;
                0 6length*Iz 0 0 0 4length^2*Iz 0 -6length*Iz 0 0 0 2length^2*Iz;
                -A*length^2 0 0 0 0 0 A*length^2 0 0 0 0 0;
                0 -12Iz 0 0 0 -6length*Iz 0 12Iz 0 0 0 -6length*Iz;
                0 0 -12Iy 0 6length*Iy 0 0 0 12Iy 0 6length*Iy 0;
                0 0 0 -G*J*length^2/E 0 0 0 0 0 G*J*length^2/E 0 0;
                0 0 -6length*Iy 0 2length^2*Iy 0 0 0 6length*Iy 0 4length^2*Iy 0;
                0 6length*Iz 0 0 0 2length^2*Iz 0 -6length*Iz 0 0 0 4length^2*Iz
                ]

            element.k = Symmetric(R' * k * R)
        elseif dims == 2
            l, m = (element.posEnd .- element.posStart) ./ element.length

            # rotation matrix
            R = [l m 0 0 0 0;
                -m l 0 0 0 0;
                0 0 1 0 0 0;
                0 0 0 l m 0;
                0 0 0 -m l 0;
                0 0 0 0 0 1]
            element.R = R

            # stiffness matrix
            k = zeros(6,6)

            p = [element.A / element.length, 0, 0, -element.A/element.length, 0, 0]
            v = [0, 12 * element.Iz/element.length^3, 6 * element.Iz/element.length^2, 0, -12element.Iz/element.length^3, 6element.Iz/element.length^2]
            m = [0, 6element.Iz/element.length^2, 4element.Iz/element.length, 0, -6element.Iz/element.length^2, 2element.Iz/element.length]

            k[1,:] = p
            k[2,:] = v
            k[3,:] = m
            k[4,:] = -p
            k[5,:] = -v
            k[6,:] = [0, 6element.Iz/element.length^2, 2element.Iz/element.length, 0, -6element.Iz/element.length^2, 4element.Iz/element.length]

            k .*= element.E

            element.k = Symmetric(R' * k * R)
        else
            error("dims must be 2 or 3")
        end
    else
        error("type must be :truss or :frame")
    end

    if return_k
        return element.k
    end
end

function K(elements::Array{Element}, nDOFS::Int)
    K = spzeros(nDOFS, nDOFS)

    for element in elements
        K[element.dofIndex, element.dofIndex] .+= element.k
    end

    return K
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

function postProcess(structure::Structure; scaleFactor = :auto)

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

function reactions(structure::Structure)
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


function addNodeElements(elements::Vector{Element}, nodes::Vector{Node})
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

function addNodeLoads(loads::Vector{Load}, nodes::Vector{Node})
   for i = 1:length(loads)
    nodes[loads[i].index].load = loads[i].load 
   end
end

function analyze(structure::Structure; forceK = false, SF = :auto)
    if !isdefined(structure, :F)
        error("No loads defined.")
    end

    # associated elements to each node
    addNodeElements(structure.elements, structure.nodes)
    # associated load vectors to individual nodes
    addNodeLoads(structure.loads, structure.nodes)

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
    nodeGlobalIndex(structure.nodes)

    #create elemental DOF indices
    dofExpander(structure.elements, structure.nodes)

    #create elemental stiffness matrices
    k_elemental.(structure.elements, structure.dims)

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
    reactions(structure)

    #post processing
    postProcess(structure; scaleFactor = SF)
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



```
Typical DOF settings:
dims : 2 or 3
type : :truss or :frame
fixity : :x/y/zfree, x/y/zfixed, :free, :fixed, :pinfixed
```
function dofMaker(dims::Int64, type::Symbol, fixity::Symbol)
    if dims == 2
        if type == :truss
            if fixity == :free
                return [true, true]
            elseif fixity == :fixed
                return [false, false]
            elseif fixity == :xfree
                return [true, false]
            elseif fixity == :yfree
                return [false, true]
            else
                error("unknown fixity, use :free, :fixed, :xfree, :yfree")
            end
        elseif type == :frame
            if fixity == :free
                return [true, true, true]
            elseif fixity == :pinfixed
                return [false, false, true]
            elseif fixity == :fixed
                return [false, false, false]
            elseif fixity == :xfree
                return [true, false, true]
            elseif fixity == :yfree
                return [false, true, false]
            else
                error("unknown fixity, use :free, :pinfixed, :fixed, :xfree, :yfree")
            end
        else
            error("type must be :truss or :frame")
        end
    elseif dims == 3
        if type == :truss
            if fixity == :free
                return [true, true, true]
            elseif fixity == :fixed
                return [false, false, false]
            elseif fixity == :zfixed
                return [true, true, false]
            elseif fixity == :xfixed
                return [false, true, true]
            elseif fixity == :yfixed
                return [true, false, true]
            elseif fixity == :xfree
                return [true, false, false]
            elseif fixity == :yfree
                return [false, true, false]
            elseif fixity == :zfree
                return [false, false, true]
            else
                error("unknown fixity, use :free, :fixed, :x/y/zfixed, :x/y/zfree")
            end
        elseif type == :frame
            if fixity == :free
                return [true, true, true, true, true, true]
            elseif fixity == :pinfixed
                return [false, false, false, true, true, true]
            elseif fixity == :fixed
                return [false, false, false, false, false, false]
            elseif fixity == :zfixed
                return [true, true, false, true, true, true]
            elseif fixity == :xfixed
                return  [false, true, true, true, true, true]
            elseif fixity == :yfixed
                return [true, false, true, true, true, true]
            elseif fixity == :xfree
                return [true, false, false, true, true, true]
            elseif fixity == :yfree
                return [false, true, false, true, true, true]
            elseif fixity == :zfree
                return [false, false, true, true, true, true]
            else
                error("unknown fixity, use :free, :pinfixed, :fixed, :x/y/zfixed, :x/y/zfree")
            end
        else
            error("type must be :truss or :frame")
        end
    else
        error("dims must be 2 or 3")
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

    analyze(structure)

    iter = 0
    while iter < maxIters
        for element in structure.elements
            element.A = trussSizer(element.axialForce, k, element.length, element.E, dDratio, maxStress, minArea)
        end
        analyze(structure; forceK = true)
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
    analyze(newstructure)

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
        analyze(newstructure; forceK = true)

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