# Rebuilds the DiffAnalysis_2024 publication structures from their serialized
# definitions in fixtures_diffanalysis.jl using ONLY plain Asap — no
# DiffAnalysis/AsapOptim/Makie needed. Used by characterization.jl to verify
# that current Asap reproduces the published forward-analysis results.

using Asap

"""
    rebuild_diffanalysis_model(def::Dict) -> TrussModel or Model

Reconstruct a publication structure from its serialized definition:
nodes (position + DOF fixity), elements (connectivity, section, roll angle Ψ,
release), and loads. The returned model is unsolved.
"""
function rebuild_diffanalysis_model(def::Dict)
    istruss = def["kind"] == "truss"

    nodes = map(def["nodes"]) do n
        pos = Vector{Float64}(n["pos"])
        dof = Vector{Bool}(n["dof"])
        istruss ? TrussNode(pos, dof) : Node(pos, dof)
    end

    elements = map(def["elements"]) do e
        s = Vector{Float64}(e["section"])
        if istruss
            TrussElement(nodes[e["i"]], nodes[e["j"]], TrussSection(s[1], s[2], s[3]))
        else
            sec = Section(s[1], s[2], s[3], s[4], s[5], s[6], s[7])
            el = Element(nodes[e["i"]], nodes[e["j"]], sec; release=Symbol(e["release"]))
            el.Ψ = e["psi"]
            el
        end
    end

    loads = map(def["loads"]) do L
        v = haskey(L, "value") ? Vector{Float64}(L["value"]) : Float64[]
        if L["type"] == "nodeforce"
            NodeForce(nodes[L["i"]], v)
        elseif L["type"] == "nodemoment"
            NodeMoment(nodes[L["i"]], v)
        elseif L["type"] == "lineload"
            LineLoad(elements[L["i"]], v)
        elseif L["type"] == "pointload"
            PointLoad(elements[L["i"]], L["position"], v)
        else
            error("unhandled load type $(L["type"])")
        end
    end

    return istruss ? TrussModel(nodes, elements, loads) : Model(nodes, elements, loads)
end
