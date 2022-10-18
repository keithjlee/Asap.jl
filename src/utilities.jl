```
Typical DOF settings:
dims : 2 or 3
type : :truss or :frame
fixity : :x/y/zfree, x/y/zfixed, :free, :fixed, :pinfixed
```
dofDict = Dict(
    (2, :truss, :free) => [true, true],
    (2, :truss, :fixed) => [false, false],
    (2, :truss, :xfree) => [true, false],
    (2, :truss, :yfree) => [false, true],
    (2, :frame, :free) => [true, true, true],
    (2, :frame, :pinfixed) => [false, false, true],
    (2, :frame, :fixed) => [false, false, false],
    (2, :frame, :xfree) => [true, false, true],
    (2, :frame, :yfree) => [false, true, false],
    (3, :truss, :free) => [true, true, true],
    (3, :truss, :fixed) => [false, false, false],
    (3, :truss, :zfixed) => [true, true, false],
    (3, :truss, :xfixed) => [false, true, true],
    (3, :truss, :yfixed) => [true, false, true],
    (3, :truss, :xfree) => [true, false, false],
    (3, :truss, :yfree) => [false, true, false],
    (3, :truss, :zfree) => [false, false, true],
    (3, :frame, :free) => [true, true, true, true, true, true],
    (3, :frame, :pinfixed) => [false, false, false, true, true, true],
    (3, :frame, :fixed) => [false, false, false, false, false, false],
    (3, :frame, :xfixed) => [false, true, true, true, true, true],
    (3, :frame, :yfixed) => [true, false, true, true, true, true],
    (3, :frame, :zfixed) => [true, true, false, true, true, true]
)

"""
Common DOF scenarios
"""
function dofMaker(dims::Int64, type::Symbol, fixity::Symbol)
    return dofDict[(dims, type, fixity)]
end

# Custom indexing
function Base.getindex(elements::Vector{Element}, i::Symbol)
    return [element for element in elements if element.id == i]
end

function Base.findall(elements::Vector{Element}, i::Symbol)
    return findall([element.id == i for element in elements])
end

function Base.getindex(nodes::Vector{Node}, i::Symbol)
    return [node for node in nodes if node.id == i]
end

function Base.findall(nodes::Vector{Node}, i::Symbol)
    return findall([node.id == i for node in nodes])
end

function Base.getindex(loads::Vector{Load}, i::Symbol)
    return [load for load in loads if load.id == i]
end

function Base.findall(loads::Vector{Load}, i::Symbol)
    return findall([load.id == i for load in loads])
end

function Base.getindex(structure::Structure, i::Symbol)
    return structure.nodes[i], structure.elements[i], structure.loads[i]
end

"""
convert Karamba model to Asap structure.
HIGHLY LIKELY TO FAIL!!!
"""
function karamba2asap(filename::String)
    if filename[end-4:end] !== ".json"
        error("file must be .json")
    end

    file = JSON.parsefile(filename)


    supportIDs = Vector{Int}()
    supportConditions = Vector{Vector{Bool}}()

    for sup in file["Supports"]["\$values"]
        push!(supportIDs, sup["NodeInd"] .+ 1)

        dofs = sup["Condition"]["\$values"]
        push!(supportConditions, [!Bool(x) for x in dofs])
    end
    supportDict = Dict(zip(supportIDs, supportConditions))

    ndof_node = length(supportConditions[1])

    nNodes = length(file["Nodes"]["\$values"])
    nodes = Vector{Node}()
    for i = 1:nNodes
        node = file["Nodes"]["\$values"][i]
        x = node["Pos"]["X"]
        y = node["Pos"]["Y"]
        z = node["Pos"]["Z"]

        if any(i .== supportIDs)
            fixity = supportDict[i]
        else
            fixity = repeat([true], ndof_node)
        end


        push!(nodes, Node([x, y, z], fixity))
        
    end

    nElements = length(file["Elems"]["\$values"])
    elements = Vector{Element}()

    karambaElements = file["Elems"]["\$values"]
    croIds = []

    for elem in karambaElements
        try 
            push!(croIds, elem["BuilderElement"]["CroSec"]["\$id"])
        catch
            push!(croIds, nothing)
        end
    end

    materialIDs = []
    for elem in karambaElements
        try 
            push!(materialIDs, elem["BuilderElement"]["CroSec"]["Material"]["\$id"])
        catch
            push!(materialIDs, nothing)
        end
    end

    elements = Vector{Element}()

    for elem in karambaElements
        ids = Int.(elem["NodeInds"]["\$values"] .+ 1)

        if haskey(elem["BuilderElement"]["CroSec"], "\$ref")
            croRefID = findfirst(croIds .== elem["BuilderElement"]["CroSec"]["\$ref"])

            croRefElem = karambaElements[croRefID]
            croRefSection = croRefElem["BuilderElement"]["CroSec"]

            A = croRefSection["A"]
            Istrong = croRefSection["Iyy"]
            Iweak = croRefSection["Izz"]
            Itwist = croRefSection["Ipp"]

            if haskey(croRefSection["Material"], "\$ref")
                matRefID = findfirst(materialIDs .== croRefSection["Material"]["\$ref"])
                matRef = karambaElements[matRefID]["BuilderElement"]["CroSec"]["Material"]
                E = matRef["E"]
                G = matRef["G12"]
            else
                E = croRefSection["Material"]["E"]
                G = croRefSection["Material"]["G12"]
            end
        else
            croRefSection = elem["BuilderElement"]["CroSec"]
            A = croRefSection["A"]
            Istrong = croRefSection["Iyy"]
            Iweak = croRefSection["Izz"]
            Itwist = croRefSection["Ipp"]

            if haskey(croRefSection["Material"], "\$ref")
                matRefID = findfirst(materialIDs .== croRefSection["Material"]["\$ref"])
                matRef = karambaElements[matRefID]["BuilderElement"]["CroSec"]["Material"]
                E = matRef["E"]
                G = matRef["G12"]
            else
                E = croRefSection["Material"]["E"]
                G = croRefSection["Material"]["G12"]
            end
        end

        push!(elements, Element(nodes, ids, E, A, G, Istrong, Iweak, Itwist))
    end

    structure = Structure(nodes, elements)

    return structure
end
