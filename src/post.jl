mutable struct GeometricNode
    node::AbstractNode
    position::Vector{Float64}
    displacedPosition::Vector{Float64}

    function GeometricNode(node::AbstractNode)
        gnode = new(node, node.position)
        gnode.displacedPosition = zeros(3)

        return gnode
    end

    function GeometricNode(node::GeometricNode)
        gnode = new(node.node, node.position)
        gnode.displacedPosition = zeros(3)

        return gnode
    end
end

function displace!(node::GeometricNode, factor::Union{Float64, Int64})
    node.displacedPosition = node.node.position .+ node.node.displacement[1:3] .* factor
end

mutable struct GeometricElement
    element::AbstractElement
    nodeStart::GeometricNode
    nodeEnd::GeometricNode
    positions::Vector{Vector{Float64}}
    nsegments::Int64
    displacedPositions::Vector{Vector{Float64}}
    midpoints::Vector{Vector{Float64}}
    LCS::Vector{Vector{Float64}}


    """
    Displaced bending element
    """
    function GeometricElement(nodes::Vector{GeometricNode}, element::Element; nsegments::Int64 = 20)
        
        gnstart = nodes[element.nodeIDs[1]]
        gnend = nodes[element.nodeIDs[2]]
        undisplaced = [gnstart.position, gnend.position]

        gelement = new(element, gnstart, gnend, undisplaced, nsegments)
        gelement.midpoints = repeat([(gnstart.position .+ gnend.position) ./ 2], 3)
        gelement.LCS = element.LCS

        return gelement
    end

    """
    Truss element
    """
    function GeometricElement(nodes::Vector{GeometricNode}, element::TrussElement; nsegments::Int64 = 20)
        
        gnstart = nodes[element.nodeIDs[1]]
        gnend = nodes[element.nodeIDs[2]]
        undisplaced = [gnstart.position, gnend.position]

        gelement = new(element, gnstart, gnend, undisplaced, 1)
        gelement.midpoints = repeat([(gnstart.position .+ gnend.position) ./ 2], 3)
        gelement.LCS = element.LCS

        return gelement
    end
end

"""
try replacing all instances of .displacedPosition with regular .position
"""
function displace!(element::GeometricElement, factor::Union{Float64, Int64})
    n1 = element.nodeStart
    n2 = element.nodeEnd

    if typeof(element.element) == TrussElement
        element.displacedPositions = [n1.displacedPosition, n2.displacedPosition]
    else
        u = [n1.node.displacement; n2.node.displacement]

        posStart = n1.displacedPosition
        L = element.element.length
        LCS = element.element.LCS
        r = element.element.R

        element.displacedPositions = disp(n1.position, n2.displacedPosition, u, L, r, LCS, element.nsegments, factor)
    end

end

mutable struct GeometricLoad
    starts::Vector{Vector{Float64}}
    vecs::Vector{Vector{Float64}}
    scalefactor::Union{Int64, Float64}

    function GeometricLoad(load::NodeForce; scalefactor = 10, n = 20)
        vec = normalize(load.value) .* scalefactor
        
        return new([load.node.position], [vec], scalefactor)
    end

    function GeometricLoad(load::NodeMoment; scalefactor = 10, n = 20)
        return new(Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), scalefactor)
    end

    function GeometricLoad(load::LineLoad; scalefactor = 10, n = 20)
        vec = normalize(load.value) .* scalefactor

        increment = load.element.length / (n-2)
        starts = Vector{Vector{Float64}}()
        dirs = Vector{Vector{Float64}}()

        for i = 1:(n-2)
            step = increment * i
            push!(starts, load.element.posStart .+ step * load.element.LCS[1])
            push!(dirs, vec)
        end

        return new(starts, dirs, scalefactor)
    end

    function GeometricLoad(load::GravityLoad; scalefactor = 10, n = 20)
        vec = normalize(load.value .* [0., 0., -1.]) .* scalefactor

        increment = load.element.length / (n-2)
        starts = Vector{Vector{Float64}}()
        dirs = Vector{Vector{Float64}}()

        for i = 1:n-2
            step = increment * i
            push!(starts, load.element.posStart .+ step * load.element.LCS[1])
            push!(dirs, vec)
        end

        return new(starts, dirs, scalefactor)
    end

    function GeometricLoad(load::PointLoad; scalefactor = 10, n = 20)
        vec = normalize(load.value) .* scalefactor
        position = load.element.posStart .+ load.position .* load.element.length .* load.element.LCS[1] 

        return new([position], [vec], scalefactor)
    end
end

mutable struct Geometry
    nodes::Vector{GeometricNode}
    elements::Vector{GeometricElement}
    loads::Vector{GeometricLoad}
    dispfactor::Union{Int64, Float64}
    loadfactor::Union{Int64, Float64}

    function Geometry(model::AbstractModel)
        nodes = [GeometricNode(node) for node in model.nodes]
        elems = [GeometricElement(nodes, element) for element in model.elements]
        loads = [GeometricLoad(load) for load in model.loads]
        factor = 1

        displace!.(nodes, factor)
        displace!.(elems, factor)

        return new(nodes, elems, loads, factor, factor)
    end

    function Geometry(model::AbstractModel, dispfactor::Union{Int64, Float64}, loadfactor::Union{Int64, Float64}; n::Int64 = 20)
        nodes = [GeometricNode(node) for node in model.nodes]
        elems = [GeometricElement(nodes, element; nsegments = n) for element in model.elements]
        loads = [GeometricLoad(load; scalefactor = loadfactor, n = n) for load in model.loads]

        displace!.(nodes, dispfactor)
        displace!.(elems, dispfactor)

        return new(nodes, elems, loads, dispfactor, loadfactor)
    end
end

function updateFactor!(geo::Geometry, factor::Union{Int64, Float64})
    geo.dispfactor = factor
    displace!.(geo.nodes, factor)
    displace!.(geo.elements, factor)
end


"""
Local coordinate system of element
"""
function lcs(element::GeometricElement, Ψ; tol = 0.001)

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

    return [xvec, yvec, zvec]
end

"""
Assumed analyzed
"""
iLocaly = [2, 6, 8, 12]
iLocalz = [3, 5, 9, 11]

iLocalx1 = [1, 2, 7, 8]
iLocalx2 = [1, 3, 7, 9]
function disp(model::Model, element::Element, n::Int64, factor::Union{Int64,Float64})
    
    #x increment
    xrange = collect(range(0, element.length, length = n))
    
    #displacement vector in LCS
    ulocal = element.R * model.u[element.globalID]

    #end displacements w/r/t y displacement
    ulocaly = ulocal[iLocaly]

    #end displacements w/r/t z displacements
    ulocalz = ulocal[iLocalz]

    #axial displacements w/r/t xy displacements
    ulocalx1 = ulocal[iLocalx1]

    #axial displacements w/r/t xz displacements
    ulocalx2 = ulocal[iLocalx2]

    #shape function
    shapeFunction = vcat([N(i, element.length) for i in xrange]...)
    shapeFunctionAxial = vcat([Naxial(i, element.length) for i in xrange]...)

    #shift factors
    xrange1 = shapeFunctionAxial * ulocalx1 * factor
    xrange2 = shapeFunctionAxial * ulocalx2 * factor
    xrange = xrange1 .+ xrange2
    yrange = shapeFunction * ulocaly * factor
    zrange = shapeFunction * ulocalz * factor

    #shift postiions
    xshift = [x * element.LCS[1] for x in xrange]
    yshift = [y * element.LCS[2] for y in yrange]
    zshift = [z * element.LCS[3] for z in zrange]

    fullshift = xshift .+ yshift .+ zshift

    return [element.posStart .+ shift for shift in fullshift]
end

function disp(posStart::Vector{Float64}, u::Vector{Float64}, L::Float64, R::Matrix{Float64}, LCS::Vector{Vector{Float64}}, n::Int64, factor::Union{Int64,Float64})
    
    
    
    #displacement vector in LCS
    ulocal = R * u

    Lnew = L - (ulocal[1] - ulocal[7])
    
    #x increment
    xrange = collect(range(0, Lnew, length = n))

    #end displacements w/r/t y displacement
    ulocaly = ulocal[iLocaly]

    #end displacements w/r/t z displacements
    ulocalz = ulocal[iLocalz]

    #shape function
    shapeFunction = vcat([N(i, L) for i in xrange]...)

    #shift factors
    yrange = shapeFunction * ulocaly * factor
    zrange = shapeFunction * ulocalz * factor

    ####TEST
    # shapeFunctionAxial = [Naxial(i, L) for i in xrange]
    
    # #axial displacements w/r/t xy displacements
    # ulocalx1 = ulocal[iLocalx1]

    # #axial displacements w/r/t xz displacements
    # ulocalx2 = ulocal[iLocalx2]

    # #shift factors
    # xyresult = vcat([sf * ulocalx1 * factor for sf in shapeFunctionAxial]...)
    # xzresult = vcat([sf * ulocalx2 * factor for sf in shapeFunctionAxial]...)

    # xy_x = xyresult[1:2:end] .* L
    # xy_y = xyresult[2:2:end]

    # xz_x = xzresult[1:2:end]
    # xz_z = xzresult[2:2:end]


    # yrange .+= xy_y
    # zrange .+= xz_z
    ###TEST

    #shift postiions
    xshift = [x * LCS[1] for x in xrange]
    yshift = [y * LCS[2] for y in yrange]
    zshift = [z * LCS[3] for z in zrange]

    fullshift = xshift .+ yshift .+ zshift

    return [posStart .+ shift for shift in fullshift]
end

function disp(posStart::Vector{Float64}, posEnd::Vector{Float64}, u::Vector{Float64}, L::Float64, R::Matrix{Float64}, LCS::Vector{Vector{Float64}}, n::Int64, factor::Union{Int64,Float64})
    
 
    #displacement vector in LCS
    ulocal = R * u
    #x increment
    Lnew = L - (ulocal[1] - ulocal[7])

    xrange = collect(range(0, Lnew, length = n))
    #end displacements w/r/t y displacement
    ulocaly = ulocal[iLocaly]

    #end displacements w/r/t z displacements
    ulocalz = ulocal[iLocalz]

    #shape function
    shapeFunction = vcat([N(i, L) for i in xrange]...)

    #shift factors
    yrange = shapeFunction * ulocaly * factor
    zrange = shapeFunction * ulocalz * factor

    #shift postiions
    # xshift = [x * LCS[1] for x in xrange]
    xshift = [x * LCS[1] for x in xrange]
    yshift = [y * LCS[2] for y in yrange]
    zshift = [z * LCS[3] for z in zrange]

    fullshift = xshift .+ yshift .+ zshift


    slack = posEnd .- (posStart .+ fullshift[end])
    slacksingle = slack ./ (length(fullshift) - 1)


    # return [posStart .+ shift for shift in fullshift]
    d = [posStart .+ shift .+ slacksingle for shift in fullshift[2:end-1]]
    p = [posStart .+ fullshift[1]]

    push!(p, d...)
    push!(p, posEnd)

    return p
end

function makeR(straightline::Vector{Float64}, shift::Vector{Float64})
    ang = round(dot(shift, straightline) / norm(shift) / norm(straightline))
    t = acos(ang) #angle 
    s = normalize(cross(shift, straightline)) #axis

    scale = norm(straightline) / norm(shift)

    R = AngleAxis(t, s...) .* scale
end