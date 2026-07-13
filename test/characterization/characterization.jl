# Characterization tests: assert Asap reproduces the exact numerics pinned in
# fixtures.jl (generated from v0.2.2). This suite is the regression contract
# for the v1.0 modernization — see docs/MODERNIZATION.md.
#
# Where the modernization legitimately changes behavior (e.g. GravityLoad →
# SelfWeight), update the corresponding test deliberately and document why.

using Asap
using Test
using LinearAlgebra

include(joinpath(@__DIR__, "models.jl"))
include(joinpath(@__DIR__, "fixtures.jl"))

# Tight tolerance: same algebra should reproduce to machine-level precision.
# rtol chosen to absorb BLAS/solver ordering differences across platforms.
approx(a, b) = isapprox(a, b; rtol=1e-9, atol=1e-11)
approx(a::Vector{<:Vector}, b::Vector{<:Vector}) = length(a) == length(b) && all(approx.(a, b))

@testset "Characterization (pinned v0.2.x numerics)" begin

    solved_models = [
        "truss_2d" => char_truss_2d,
        "frame_2d" => char_frame_2d,
        "truss_3d" => char_truss_3d,
        "frame_3d" => char_frame_3d,
        "cantilever_tip" => char_cantilever_tip,
        "cantilever_midpoint" => char_cantilever_midpoint,
        ["portal_$r" => (() -> char_portal_frame(r)) for r in CHAR_RELEASES]...,
    ]

    @testset "solve: $name" for (name, builder) in solved_models
        model = builder()
        solve!(model)
        @test approx(model.u, FIXTURES["$name/u"])
        @test approx(model.reactions, FIXTURES["$name/reactions"])
        @test approx(model.compliance, FIXTURES["$name/compliance"])
        @test approx([el.forces for el in model.elements], FIXTURES["$name/element_forces"])
    end

    @testset "element matrices + FEFs: $r" for r in CHAR_RELEASES
        model = char_fef_element(r)
        process!(model)
        el = model.elements[1]

        @test approx(el.R, FIXTURES["fef_$r/R"])
        @test approx(el.LCS, FIXTURES["fef_$r/LCS"])
        @test approx(el.length, FIXTURES["fef_$r/length"])
        @test approx(Asap.local_K(el), FIXTURES["fef_$r/local_K"])
        @test approx(el.K, FIXTURES["fef_$r/global_K"])
        @test approx(el.Q, FIXTURES["fef_$r/Q"])

        for load in model.loads
            kind = load isa LineLoad ? "lineload" : "pointload"
            @test approx(Asap.q_local(load), FIXTURES["fef_$r/q_local/$kind"])
            @test approx(Asap.q(load), FIXTURES["fef_$r/q/$kind"])
        end
    end

    # The AsapToolkit InternalForces oracles (consumed in Phase 3) are
    # generated separately (generate_toolkit_fixtures.jl); just validate the
    # committed file parses and is populated.
    @testset "toolkit fixtures file" begin
        include(joinpath(@__DIR__, "fixtures_toolkit.jl"))
        @test TOOLKIT_FIXTURES isa Dict{String,Any}
        @test length(TOOLKIT_FIXTURES) >= 80
    end

    # Known bug: q_local(::GravityLoad) references nonexistent fields
    # (load.element.ρ and a bare `element`), so any model carrying a
    # GravityLoad throws during processing. When SelfWeight lands (Phase 2),
    # this flips to an unexpected pass — replace it with a real test then.
    @testset "GravityLoad (documented bug)" begin
        model = char_gravity_model()
        @test_broken (solve!(model); true)
    end
end
