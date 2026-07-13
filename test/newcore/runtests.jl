# New-core (Phase 1, dark development) test entry point.
# Included from the main test/runtests.jl; can also be run standalone:
#     julia --project=. -e 'using Test; include("test/newcore/runtests.jl")'

using Test
using Asap            # legacy module — compatibility oracles (fixDict etc.)
using StaticArrays

include("harness.jl")

@testset "New core (Phase 1)" begin
    include("test_materials.jl")
    include("test_nodes_ends.jl")
end
