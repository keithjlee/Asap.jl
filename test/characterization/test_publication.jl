# DiffAnalysis_2024 publication structures, rebuilt with the v1.0 API and
# verified against the pinned legacy results (fixtures_diffanalysis.jl,
# generated under the publication's Asap 0.2.1 environment).
#
# DOF-space mapping as in the parity tests: legacy truss models used 3
# DOFs/node (v1.0 always 6, rotations inactive); legacy truss element forces
# were 2-vectors (v1.0 uniform 12-vectors).
#
# Tolerance note: models containing off-midspan PointLoads with an axial
# component would reflect the documented FEF deviation — the publication
# structures carry only NodeForces and midspan/transverse element loads, so
# exact agreement is expected.

"""
Rebuild a serialized publication structure with the v1.0 API.
"""
function rebuild_publication_model(def::Dict)
    istruss = def["kind"] == "truss"

    nodes = map(def["nodes"]) do n
        pos = Vector{Float64}(n["pos"])
        dof = Vector{Bool}(n["dof"])
        fixity = istruss ? vcat(dof, [true, true, true]) : dof
        Node(pos, fixity)
    end

    elements = map(def["elements"]) do e
        s = Vector{Float64}(e["section"])
        if istruss
            mat = Material(s[2], 1.0, s[3], 0.3)          # E, (G unused), ρ
            TrussElement(nodes[e["i"]], nodes[e["j"]], Section(mat, s[1]))
        else
            mat = Material(s[2], s[3], s[7], 0.3)          # E, G, ρ
            sec = Section(mat, s[1], s[4], s[5], s[6])     # A, Ix, Iy, J
            el = FrameElement(nodes[e["i"]], nodes[e["j"]], sec;
                release=Symbol(e["release"]), rollangle=e["psi"])
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

    return Model(nodes, collect(AbstractElement{Float64}, elements),
        collect(AbstractLoad{Float64}, loads))
end

@testset "DiffAnalysis publication models (golden source)" begin
    include(joinpath(@__DIR__, "fixtures_diffanalysis.jl"))
    model_keys = sort!(collect(filter(k -> endswith(k, "/model"), keys(DIFFANALYSIS_FIXTURES))))
    @test !isempty(model_keys)

    @testset "$key" for key in model_keys
        def = DIFFANALYSIS_FIXTURES[key]
        istruss = def["kind"] == "truss"
        ndl = istruss ? 3 : 6

        model = rebuild_publication_model(def)
        solve!(model)
        res = model.results
        nn = length(model.nodes)

        @test res.u ≈ expand_legacy(Vector{Float64}(def["u"]), ndl, nn) rtol = 1e-8 atol = 1e-10
        @test res.reactions ≈ expand_legacy(Vector{Float64}(def["reactions"]), ndl, nn) rtol = 1e-8 atol = 1e-8
        @test res.compliance ≈ def["compliance"] rtol = 1e-8

        legacy_forces = def["element_forces"]
        for (i, el) in enumerate(model.elements)
            @test element_forces(res, el) ≈
                  expand_legacy_forces(Vector{Float64}(legacy_forces[i])) rtol = 1e-7 atol = 1e-7
        end
    end
end
