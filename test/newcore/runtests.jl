# New-core (Phase 1, dark development) test entry point.
# Included from the main test/runtests.jl; can also be run standalone:
#     julia --project=. -e 'using Test; include("test/newcore/runtests.jl")'

using Test
using Asap            # legacy module — compatibility oracles (fixDict etc.)
using StaticArrays
using LinearAlgebra

include("harness.jl")
include(joinpath(@__DIR__, "..", "characterization", "fixtures.jl"))  # pinned oracles

const CHAR_RELEASES = [:fixedfixed, :fixedfree, :freefixed, :freefree, :joist]

@testset "New core (Phase 1)" begin
    include("test_materials.jl")
    include("test_nodes_ends.jl")
    include("test_kernels.jl")
    include("test_elements.jl")
    include("test_pipeline.jl")
    include("test_variable.jl")
    include("test_ad.jl")
end
