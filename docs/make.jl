# Build the Asap.jl documentation site.
#
# Local build:
#   julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
#   julia --project=docs docs/make.jl
# then open docs/build/index.html (prettyurls are off outside CI, so direct
# file browsing works).

using Documenter, Asap
using Zygote, LinearSolve   # activate AsapChainRulesExt / AsapLinearSolveExt for guide @examples

makedocs(
    sitename = "Asap.jl",
    authors = "Keith J. Lee",
    modules = [Asap],
    checkdocs = :exports,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://keithjlee.github.io/Asap.jl",
        edit_link = "main",
        size_threshold = 400 * 2^10,
        size_threshold_warn = 250 * 2^10,
    ),
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "Materials and sections" => "guide/materials-sections.md",
            "Nodes and supports" => "guide/nodes-supports.md",
            "Elements" => "guide/elements.md",
            "Loads and springs" => "guide/loads-springs.md",
            "Analysis and results" => "guide/analysis-results.md",
            "Internal forces and deflections" => "guide/internal-forces.md",
            "Load cases and combinations" => "guide/load-cases.md",
            "Differentiable analysis" => "guide/differentiable.md",
            "Force density method" => "guide/fdm.md",
            "Structure generators" => "guide/generators.md",
            "Geometry extraction" => "guide/geometry.md",
            "Solver backends" => "guide/solvers.md",
        ],
        "Migrating from v0.2" => "migration.md",
        "API reference" => [
            "Materials and sections" => "api/materials.md",
            "Nodes" => "api/nodes.md",
            "Elements" => "api/elements.md",
            "Loads and springs" => "api/loads.md",
            "Model" => "api/model.md",
            "Analysis and solving" => "api/analysis.md",
            "Results and recovery" => "api/results.md",
            "Differentiable path" => "api/differentiable.md",
            "Force density method" => "api/fdm.md",
            "Generators" => "api/generators.md",
            "Geometry" => "api/geometry.md",
        ],
    ],
)

deploydocs(
    repo = "github.com/keithjlee/Asap.jl",
    devbranch = "main",
    push_preview = false,
)
