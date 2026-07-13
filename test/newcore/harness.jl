# Dark-development harness for the Phase 1 new core (docs/MODERNIZATION.md).
#
# The new-core source files live at their FINAL src/ paths but are not yet
# included by the Asap module — they are wrapped here in a scratch module so
# they can be developed and tested against the Phase 0 oracles without
# touching the legacy code. When the full pipeline passes characterization,
# src/Asap.jl flips to these includes and this harness dissolves.

module AsapNext

using LinearAlgebra, SparseArrays, StaticArrays

const SRC = joinpath(@__DIR__, "..", "..", "src")

include(joinpath(SRC, "Materials", "materials.jl"))
include(joinpath(SRC, "Materials", "sections.jl"))
include(joinpath(SRC, "ShowMethods.jl"))

end # module AsapNext
