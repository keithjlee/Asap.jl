# Captures golden reference data from the DiffAnalysis_2024 publication repo
# (Lee, Huang & Mueller — "A differentiable structural analysis framework for
# high-performance design optimization"), the cleaned-up Asap + AsapOptim-style
# pipeline in action. Run with the Julia version matching the publication
# Manifest (julia_version = 1.11.x):
#
#     julia +1.11 test/characterization/generate_diffanalysis_fixtures.jl
#
# Runs each paper problem's setup script VERBATIM under the publication's own
# pinned environment (Asap 0.2.1), except for patched-out interactive lines
# (window display, figure/JSON export, an `import AsapToolkit` absent from the
# pub env) — the structural math is untouched. Captures:
#
#   1. Full model definitions (nodes/DOFs/elements/sections/releases/loads),
#      so tests can REBUILD each structure with plain Asap and compare — no
#      DiffAnalysis/Makie needed at test time.
#   2. Forward results: u, reactions, compliance, element end forces.
#   3. Differentiable-path data at the initial design point x_init: objective
#      value + Zygote gradient, constraint vector + jacobian (full if small,
#      deterministic probe products J·v / Jᵀ·w if large). These are the golden
#      AD reference for the Phase 5a AsapOptim rework.
#
# Full optimization runs are deliberately NOT pinned: endpoints depend on
# NLopt/solver versions and wall-clock termination, so they make poor fixtures.

using Pkg

const PUBDIR = joinpath(
    homedir(),
    "Library/CloudStorage/OneDrive-MassachusettsInstituteofTechnology",
    "Publications/DiffAnalysis_Structures_2024/DiffAnalysis_2024_CODE",
)
const SCRIPTS = joinpath(PUBDIR, "paper-scripts")

Pkg.activate(PUBDIR)
Pkg.instantiate()

using LinearAlgebra
using Asap

const PATCHDIR = mktempdir()

# Lines that only make sense interactively (or need packages absent from the
# pub env). Everything else — including all structural math — runs verbatim.
const DROP_PATTERNS = [
    r"^\s*import AsapToolkit",
    r"\batk\.",             # GHsave exports
    r"^\s*display\(",       # window display
    r"^\s*save\(",          # figure export into the pub repo
]

"""
Include `path` into `mod` with interactive-only lines commented out and
`@__DIR__` rewritten to the file's true directory (so data files resolve even
though the patched copy lives in a temp dir). cd's to the original directory
for the include so relative reads/writes stay inside the pub repo layout.
"""
function include_patched(mod::Module, path::String)
    dir = dirname(abspath(path))
    lines = split(read(path, String), '\n')
    patched = map(lines) do line
        any(p -> occursin(p, line), DROP_PATTERNS) ? "# [patched out] " * line : line
    end
    src = replace(join(patched, '\n'), "@__DIR__" => repr(dir))
    file = joinpath(PATCHDIR, string(hash(path), "_", basename(path)))
    write(file, src)
    cd(dir) do
        Base.include(mod, file)
    end
end

function problem_module()
    mod = Module(gensym(:pub))
    Core.eval(mod, :(using DiffAnalysis, LinearAlgebra))
    return mod
end

release_string(::Asap.Element{R}) where {R} = lowercase(String(nameof(R)))

"Serialize a solved model: full definition (rebuildable with plain Asap) + results."
function serialize_model(model)
    istruss = model isa TrussModel

    nodes = [Dict{String,Any}("pos" => copy(n.position), "dof" => collect(n.dof))
             for n in model.nodes]

    elements = map(model.elements) do e
        d = Dict{String,Any}("i" => e.nodeStart.nodeID, "j" => e.nodeEnd.nodeID)
        s = e.section
        if e isa TrussElement
            d["section"] = [s.A, s.E, s.ρ]
        else
            d["section"] = [s.A, s.E, s.G, s.Ix, s.Iy, s.J, s.ρ]
            d["psi"] = e.Ψ
            d["release"] = release_string(e)
        end
        d
    end

    loads = map(model.loads) do L
        if L isa NodeForce
            Dict{String,Any}("type" => "nodeforce", "i" => L.node.nodeID, "value" => copy(L.value))
        elseif L isa NodeMoment
            Dict{String,Any}("type" => "nodemoment", "i" => L.node.nodeID, "value" => copy(L.value))
        elseif L isa LineLoad
            Dict{String,Any}("type" => "lineload", "i" => L.element.elementID, "value" => copy(L.value))
        elseif L isa PointLoad
            Dict{String,Any}("type" => "pointload", "i" => L.element.elementID,
                "position" => L.position, "value" => copy(L.value))
        else
            error("unhandled load type $(typeof(L))")
        end
    end

    Dict{String,Any}(
        "kind" => istruss ? "truss" : "frame",
        "nodes" => nodes,
        "elements" => elements,
        "loads" => loads,
        # results
        "u" => copy(model.u),
        "reactions" => copy(model.reactions),
        "compliance" => model.compliance,
        "element_forces" => [copy(e.forces) for e in model.elements],
    )
end

# Deterministic probe vectors for large jacobians (no RNG — reproducible).
probe(n::Int) = normalize!(cos.(1.0:n) .+ 0.1)

"Pin a jacobian: full if small, J·v and Jᵀ·w probe products otherwise."
function pin_jacobian!(fixtures, key, J::AbstractMatrix)
    if length(J) <= 60_000
        fixtures["$key/full"] = Matrix(J)
    else
        nc, nx = size(J)
        fixtures["$key/Jv"] = Matrix(J) * probe(nx)
        fixtures["$key/Jtw"] = Matrix(J)' * probe(nc)
        fixtures["$key/size"] = [nc, nx]
    end
end

fixtures = Dict{String,Any}()
skipped = String[]

# ─────────────────────────────────────────────────────────────────────────────
# S4.1 — minimum-volume Warren truss (2D truss; area + spatial variables)
# ─────────────────────────────────────────────────────────────────────────────
try
    mod = problem_module()
    include_patched(mod, joinpath(SCRIPTS, "S4.1_minvolume_warren", "init_problem.jl"))
    Core.eval(mod, :(Asap.solve!(model)))
    fixtures["S4.1_warren/model"] = serialize_model(mod.model)

    # objective/constraints exactly as in gb_mma.jl
    o, doo, c, dc, x0 = Core.eval(mod, quote
        OBJ = x -> begin
            geo = GeometricProperties(x, params)
            dot(geo.A, geo.L)
        end
        CSTR = x -> begin
            res = solve_truss(x, params)
            [(-res.U[2:3:end] .- dmax); (axial_stress(res, params) .- fy)]
        end
        og = Zygote.withgradient(OBJ, x_init)
        cj = Zygote.withjacobian(CSTR, x_init)
        (og.val, og.grad[1], cj.val, cj.grad[1], copy(x_init))
    end)
    fixtures["S4.1_warren/x_init"] = x0
    fixtures["S4.1_warren/objective"] = o
    fixtures["S4.1_warren/objective_grad"] = doo
    fixtures["S4.1_warren/constraints"] = c
    pin_jacobian!(fixtures, "S4.1_warren/constraint_jac", dc)
    println("captured S4.1_warren")
catch err
    push!(skipped, "S4.1_warren: " * first(split(sprint(showerror, err), '\n')))
end

# ─────────────────────────────────────────────────────────────────────────────
# S4.2 — minimum-volume spaceframe (3D truss; area variables)
# ─────────────────────────────────────────────────────────────────────────────
try
    mod = problem_module()
    include_patched(mod, joinpath(SCRIPTS, "S4.2_minvolume_spaceframe", "init_problem.jl"))
    Core.eval(mod, :(Asap.solve!(model)))
    fixtures["S4.2_spaceframe/model"] = serialize_model(mod.model)

    o, doo, c, dc, x0 = Core.eval(mod, quote
        OBJ = x -> begin
            geo = GeometricProperties(x, params)
            dot(geo.L, geo.A)
        end
        CSTR = x -> begin
            res = solve_truss(x, params)
            [(-res.U[3:3:end] .- dmax); (axial_stress(res, params)[i_stressed_elements] .- fy)]
        end
        og = Zygote.withgradient(OBJ, x_init)
        cj = Zygote.withjacobian(CSTR, x_init)
        (og.val, og.grad[1], cj.val, cj.grad[1], copy(x_init))
    end)
    fixtures["S4.2_spaceframe/x_init"] = x0
    fixtures["S4.2_spaceframe/objective"] = o
    fixtures["S4.2_spaceframe/objective_grad"] = doo
    fixtures["S4.2_spaceframe/constraints"] = c
    pin_jacobian!(fixtures, "S4.2_spaceframe/constraint_jac", dc)
    println("captured S4.2_spaceframe")
catch err
    push!(skipped, "S4.2_spaceframe: " * first(split(sprint(showerror, err), '\n')))
end

# ─────────────────────────────────────────────────────────────────────────────
# S4.3 — freeform frame (3D FRAME model; forward analysis)
# ─────────────────────────────────────────────────────────────────────────────
try
    mod = problem_module()
    include_patched(mod, joinpath(SCRIPTS, "S4.3_freeform_frame", "initialize.jl"))
    Core.eval(mod, :(Asap.solve!(model)))
    fixtures["S4.3_frame/model"] = serialize_model(mod.model)
    println("captured S4.3_frame")
catch err
    push!(skipped, "S4.3_frame: " * first(split(sprint(showerror, err), '\n')))
end

# ─────────────────────────────────────────────────────────────────────────────
# S4.4 — embodied carbon bridge (large truss; forward analysis)
# ─────────────────────────────────────────────────────────────────────────────
try
    mod = problem_module()
    include_patched(mod, joinpath(SCRIPTS, "S4.4_embodiedcarbon_bridge", "init.jl"))
    Core.eval(mod, :(Asap.solve!(model)))
    fixtures["S4.4_bridge/model"] = serialize_model(mod.model)
    println("captured S4.4_bridge")
catch err
    push!(skipped, "S4.4_bridge: " * first(split(sprint(showerror, err), '\n')))
end

# ─────────────────────────────────────────────────────────────────────────────
# Write
# ─────────────────────────────────────────────────────────────────────────────
asap_ver = string(pkgversion(Asap))
path = joinpath(@__DIR__, "fixtures_diffanalysis.jl")
open(path, "w") do io
    println(io, "# AUTO-GENERATED by generate_diffanalysis_fixtures.jl.")
    println(io, "# Golden reference from the DiffAnalysis_2024 publication repo,")
    println(io, "# run under its pinned environment: Asap v", asap_ver,
        ", Julia v", VERSION, ".")
    println(io, "# Model definitions are rebuildable with plain Asap (see")
    println(io, "# diffanalysis_models.jl); gradient data is the Phase 5a AD oracle.")
    for s in skipped
        println(io, "# SKIPPED: ", s)
    end
    println(io, "DIFFANALYSIS_FIXTURES = Dict{String,Any}(")
    for key in sort!(collect(keys(fixtures)))
        println(io, "    ", repr(key), " => ", repr(fixtures[key]), ",")
    end
    println(io, ")")
end

println("Wrote $(length(fixtures)) fixtures to $path")
isempty(skipped) || println("Skipped:\n  " * join(skipped, "\n  "))
