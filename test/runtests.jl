# Asap v1.0 test suite.
#
# The regression backbone is the Phase 0 characterization data (exact
# numerics pinned from the legacy v0.2.x implementation and from the
# DiffAnalysis_2024 publication environment) — see docs/MODERNIZATION.md.
# Deliberate deviations from legacy behavior are documented there and
# asserted explicitly where they occur.

using Asap
using Test
using LinearAlgebra
using SparseArrays
using StaticArrays
using Zygote            # loads ChainRulesCore → activates AsapChainRulesExt
using FiniteDifferences

include("newcore/harness.jl")                               # AsapNext = Asap alias
include(joinpath(@__DIR__, "characterization", "fixtures.jl"))

const CHAR_RELEASES = [:fixedfixed, :fixedfree, :freefixed, :freefree, :joist]

@testset "Asap v1.0" begin
    include("newcore/test_materials.jl")
    include("newcore/test_nodes_ends.jl")
    include("newcore/test_kernels.jl")
    include("newcore/test_elements.jl")
    include("newcore/test_pipeline.jl")
    include("newcore/test_loads.jl")
    include("newcore/test_variable.jl")
    include("newcore/test_ad.jl")
    include("characterization/test_publication.jl")

    # AsapToolkit oracles: consumed by the Phase 3 recovery tests
    include(joinpath(@__DIR__, "characterization", "fixtures_toolkit.jl"))
    @testset "toolkit fixtures file" begin
        @test TOOLKIT_FIXTURES isa Dict{String,Any}
        @test length(TOOLKIT_FIXTURES) >= 80
    end
    include("newcore/test_recovery.jl")
    include("newcore/test_cases.jl")
    include("newcore/test_generation.jl")
    include("newcore/test_solvers.jl")
    include("newcore/test_fdm.jl")
end
