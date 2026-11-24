#=
Walkthrough A2 — a 3D space truss (double-layer grid)

Shift+enter through the blocks. Units: kN, m.

A 12 m × 9 m double-layer spaceframe roof from the `SpaceFrame` generator:
member force distribution, reactions, self-weight, and how the support
layout changes the structure's behaviour.
=#

begin
    using Asap, LinearAlgebra
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

begin
    steel = Material(200e6, 77e6, 7.85, 0.3)
    tube = Section(steel, 1.2e-3)                       # 12 cm² tube
    sf = SpaceFrame(8, 1.5, 6, 1.5, 1.2, tube; support = :corner, load = [0.0, 0.0, -4.0])
    roof = sf.model
    roof
end

#=
## Add self-weight and re-solve

Generators return a solved model; it is an ordinary `Model`, so you can add
loads and solve again. A `TrussElement` carries axial force only, so it
cannot take an element load (`SelfWeight` is for frames and shells) — lump
each member's weight ρA·g·L half to each end node instead, the standard
truss idealization.
=#
begin
    function lumped_selfweight(model; g = 9.81)
        acc = Dict{Node{Float64},Float64}()
        for el in model.elements
            W = ρA(el.section) * g * length(el)
            acc[el.nodeStart] = get(acc, el.nodeStart, 0.0) + W / 2
            acc[el.nodeEnd] = get(acc, el.nodeEnd, 0.0) + W / 2
        end
        return AbstractLoad{Float64}[NodeForce(n, [0.0, 0.0, -W]) for (n, W) in acc]
    end
    append!(roof.loads, lumped_selfweight(roof))
    solve!(roof)
    res = roof.results
    N = [axial_force(res, el) for el in roof.elements]
    # the generator loads every unsupported node (both layers), and Asap's `reaction`
    # leaves out loads applied directly at supports: compare with the load elsewhere
    supported = Set(roof.nodes[sf.isupport])
    applied = sum(l -> -l.value[3], filter(l -> !(l.node in supported), roof.loads))
    Rz = sum(i -> reaction(res, roof.nodes[i])[3], sf.isupport)
    (max_tension = maximum(N), max_compression = minimum(N), applied = applied, reactions = Rz,
     max_deflection_mm = 1000 * maximum(abs, res.u[3:6:end]))
end

begin
    f = fig_3d(roof; scale = 50.0, linecolor = :axial, title = "corner-supported spaceframe: axial force, deflection ×50")
    savefig("A2_spaceframe.png", f)
end

#=
## Support layouts

`support = :corner` vs `:x` (the two short edges at x = 0 and 12 m, so the
roof spans 12 m in X) vs `:y` (the two long edges, spanning 9 m in Y) vs
`:xy` (all edges): the same grid becomes a two-way plate on corner
supports, a one-way spanning roof, or a four-edge-supported plate.
=#
begin
    layouts = [:corner, :x, :y, :xy]
    defl = map(layouts) do s
        m = SpaceFrame(8, 1.5, 6, 1.5, 1.2, tube; support = s, load = [0.0, 0.0, -4.0]).model
        1000 * maximum(abs, m.results.u[3:6:end])
    end
    f = fig_bars(string.(layouts), defl; ylabel = "max deflection [mm]", title = "support layout")
    savefig("A2_supports.png", f)
    defl
end

#=
## Numbers for the report

`reports/A2_spaceframe_3d.pdf` prints these (see `figures.jl`, `report_data`).
=#
begin
    selfweight = sum(el -> ρA(el.section) * 9.81 * length(el), roof.elements)
    groups = [:top, :bottom, :web]
    supnodes = Set(roof.nodes[sf.isupport])
    n_loaded = count(n -> !(n in supnodes), roof.nodes)      # the generator loads every unsupported node, both layers
    sw_at_supports = sum(l -> -l.value[3], filter(l -> l.node in supnodes, lumped_selfweight(roof)); init = 0.0)
    report_data("A2";
        nx = sf.nx, ny = sf.ny, dx = sf.dx, dy = sf.dy, dz = sf.dz,
        plan_x = sf.nx * sf.dx, plan_y = sf.ny * sf.dy,
        E = steel.E, rho = steel.ρ, area = tube.A,
        n_nodes = length(roof.nodes), n_members = length(roof.elements),
        n_top = length(sf.itop), n_supports = length(sf.isupport), node_load = 4.0,
        n_loaded = n_loaded, node_load_total = 4.0 * n_loaded, selfweight_total = selfweight,
        reactions = Rz, selfweight_at_supports = sw_at_supports,
        reactions_expected = 4.0 * n_loaded + selfweight - sw_at_supports,
        max_tension = maximum(N), max_compression = minimum(N),
        max_deflection_mm = 1000 * maximum(abs, res.u[3:6:end]),
        forces = [(group = string(g), count = count(el -> el.id == g, roof.elements),
                   min = minimum(N[i] for (i, el) in enumerate(roof.elements) if el.id == g),
                   max = maximum(N[i] for (i, el) in enumerate(roof.elements) if el.id == g)) for g in groups],
        layouts = (names = string.(layouts), deflection_mm = defl,
                   supports = [length(SpaceFrame(8, 1.5, 6, 1.5, 1.2, tube; support = s, load = [0.0, 0.0, -4.0]).isupport) for s in layouts]))
end
