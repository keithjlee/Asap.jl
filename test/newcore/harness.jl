# Dark-development harness for the Phase 1 new core (docs/MODERNIZATION.md).
#
# The new-core source files live at their FINAL src/ paths but are not yet
# included by the Asap module — they are wrapped here in a scratch module so
# they can be developed and tested against the Phase 0 oracles without
# touching the legacy code. When the full pipeline passes characterization,
# src/Asap.jl flips to these includes and this harness dissolves.

module AsapNext

using LinearAlgebra, SparseArrays, StaticArrays
import ChainRulesCore

const SRC = joinpath(@__DIR__, "..", "..", "src")

include(joinpath(SRC, "Materials", "materials.jl"))
include(joinpath(SRC, "Materials", "sections.jl"))
include(joinpath(SRC, "Nodes", "node.jl"))
include(joinpath(SRC, "Elements", "end_conditions.jl"))
include(joinpath(SRC, "Elements", "kernels", "transformation.jl"))
include(joinpath(SRC, "Elements", "kernels", "stiffness.jl"))
include(joinpath(SRC, "Elements", "interface.jl"))
include(joinpath(SRC, "Elements", "frame.jl"))
include(joinpath(SRC, "Springs", "nodal_springs.jl"))
include(joinpath(SRC, "Loads", "loads.jl"))
include(joinpath(SRC, "Elements", "kernels", "fixed_end_forces.jl"))
include(joinpath(SRC, "Model", "model.jl"))
include(joinpath(SRC, "Analysis", "dofs.jl"))
include(joinpath(SRC, "Analysis", "symbolic.jl"))
include(joinpath(SRC, "Analysis", "assemble.jl"))
include(joinpath(SRC, "Results", "results.jl"))
include(joinpath(SRC, "Analysis", "solve.jl"))
include(joinpath(SRC, "Analysis", "functional.jl"))
include(joinpath(SRC, "ShowMethods.jl"))

# the AD rules (a package extension after the module flip; included directly
# during dark development — HOST resolution handles both)
include(joinpath(SRC, "..", "ext", "AsapChainRulesExt.jl"))

end # module AsapNext
