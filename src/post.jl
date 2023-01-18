mutable struct GeometricNode
    position::Vector{Float64}
    initPosition::Vector{Float64}
    displacement::Vector{Float64}
    id::Union{Symbol, Nothing}

    function GeometricNode(node::AbstractNode)
        gnode = new(node.position)
        gnode.initPosition = node.position
        gnode.displacement = zeros(6)
        gnode.id = node.id
        return gnode
    end

    function GeometricNode(node::AbstractNode, factor::Union{Float64, Int64})
        position = node.position .+ node.displacement[1:3] .* factor
        gnode = new(position)
        gnode.initPosition = node.position
        gnode.displacement = node.displacement
        gnode.id = node.id

        return gnode
    end

    function GeometricNode(node::GeometricNode, factor::Union{Float64, Int64})
        position = node.position .+ node.displacement[1:3] .* factor
        gnode = new(position)
        gnode.initPosition = node.initPosition
        gnode.displacement = node.displacement
        gnode.id = node.id

        return gnode
    end

end

mutable struct GeometricElement
    posStart::Vector{Float64}
    posEnd::Vector{Float64}
    length::Float64
    LCS::Vector{Vector{Float64}}
    R::Matrix{Float64}
    positions::Vector{Vector{Float64}}
    displacement::Vector{Float64}
    id::Union{Symbol, Nothing}
    nsegments::Int64


    """
    Displaced bending element
    """
    function GeometricElement(nodes::Vector{GeometricNode}, element::Element, factor::Union{Float64, Int64}, nsegments::Int64)
        
        pstart = nodes[element.nodeIDs[1]].position
        pend = nodes[element.nodeIDs[2]].position
        l = norm(pend .- pstart)

        gelement = new(pstart, pend, l)
        gelement.LCS = lcs(gelement, element.Ψ)
        gelement.R = element.R
        gelement.nsegments = nsegments

        gelement.displacement = [nodes[element.nodeIDs[1]].displacement; nodes[element.nodeIDs[2]].displacement]

        gelement.positions = disp(gelement.posStart, gelement.displacement, gelement.length, gelement.R, gelement.LCS, nsegments, factor)

        gelement.id = element.id

        return gelement
    end

    """
    Displaced bending element
    """
    function GeometricElement(nodes::Vector{GeometricNode}, element::GeometricElement, factor::Union{Float64, Int64}, nsegments::Int64)
        
        pstart = element.posStart
        pend = element.posEnd
        l = norm(pend .- pstart)

        gelement = new(pstart, pend, l)
        gelement.LCS = element.LCS
        gelement.R = element.R
        gelement.nsegments = nsegments

        gelement.displacement = element.displacement

        gelement.positions = disp(gelement.posStart, gelement.displacement, gelement.length, gelement.R, gelement.LCS, nsegments, factor)

        gelement.id = element.id

        return gelement
    end

    """
    Bending element
    """
    function GeometricElement(nodes::Vector{GeometricNode}, element::Element)
        
        pstart = nodes[element.nodeIDs[1]].position
        pend = nodes[element.nodeIDs[2]].position
        l = norm(pend .- pstart)

        gelement = new(pstart, pend, l)
        gelement.LCS = lcs(gelement, element.Ψ)
        gelement.R = element.R
        gelement.nsegments = 1

        gelement.positions = [pstart, pend]

        gelement.id = element.id

        return gelement
    end

    """
    Truss element
    """
    function GeometricElement(nodes::Vector{GeometricNode}, element::TrussElement)
        
        pstart = nodes[element.nodeIDs[1]].position
        pend = nodes[element.nodeIDs[2]].position
        l = norm(pend .- pstart)

        gelement = new(pstart, pend, l)
        gelement.LCS = Vector{Vector{Float64}}()
        gelement.R = element.R
        gelement.nsegments = 1

        gelement.positions = [pstart, pend]

        gelement.id = element.id

        return gelement
    end
end

mutable struct Geometry
    nodes::Vector{GeometricNode}
    elements::Vector{GeometricElement}
    displacedNodes::Vector{GeometricNode}
    displacedElements::Vector{GeometricElement}
    displacementFactor::Union{Int64, Float64}

    function Geometry(model::AbstractModel)
        nodes = [GeometricNode(node) for node in model.nodes]
        elems = [GeometricElement(nodes, element) for element in model.elements]
        dnodes = Vector{GeometricNode}()
        delems = Vector{GeometricElement}()
        factor = 1

        return new(nodes, elems, dnodes, delems, factor)
    end
    
    function Geometry(model::AbstractModel, factor::Union{Int64, Float64}, n::Int64)
        nodes = [GeometricNode(node) for node in model.nodes]
        elems = [GeometricElement(nodes, element) for element in model.elements]
        dnodes = [GeometricNode(node, factor) for node in model.nodes]
        delems = Vector{GeometricElement}()
        for element in model.elements
            if typeof(element) == Element
                push!(delems, GeometricElement(dnodes, element, factor, n))
            else
                push!(delems, GeometricElement(dnodes, element))
            end
        end

        return new(nodes, elems, dnodes, delems, factor)
    end
end

function updateFactor!(geo::Geometry, factor::Union{Int64, Float64})
    geo.displacementFactor = factor
    newdnodes = [GeometricNode(node, factor) for node in geo.nodes]
    newdelems = Vector{GeometricElement}()
    for element in geo.displacedElements
        if length(element.positions) > 2
            push!(newdelems, GeometricElement(newdnodes, element, factor, element.nsegments))
        else
            push!(newdelems, GeometricElement(newdnodes, element))
        end
    end

    geo.displacedNodes = newdnodes
    geo.displacedElements = newdelems
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
function disp(model::Model, element::Element, n::Int64, factor::Union{Int64,Float64})
    
    #x increment
    xrange = collect(range(0, element.length, length = n))
    
    #displacement vector in LCS
    ulocal = element.R * model.u[element.globalID]

    #end displacements w/r/t y displacement
    ulocaly = ulocal[iLocaly]

    #end displacements w/r/t z displacements
    ulocalz = ulocal[iLocalz]

    #shape function
    shapeFunction = vcat([N(i, element.length) for i in xrange]...)

    #shift factors
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
    
    #x increment
    xrange = collect(range(0, L, length = n))
    
    #displacement vector in LCS
    ulocal = R * u

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
    xshift = [x * LCS[1] for x in xrange]
    yshift = [y * LCS[2] for y in yrange]
    zshift = [z * LCS[3] for z in zrange]

    fullshift = xshift .+ yshift .+ zshift

    return [posStart .+ shift for shift in fullshift]
end
