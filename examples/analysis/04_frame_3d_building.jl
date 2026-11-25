#=
Walkthrough A4 — a 3D building frame: load cases, combinations, envelopes

Shift+enter through the blocks. Units: kN, m.

The `Frame` generator builds a 3D steel frame (columns, primary beams,
joists with pinned ends, braces). We replace its demo loads with real load
CASES — dead, live, wind — solve all of them against one factorization,
combine them, and envelope a beam.
=#

begin
    using Asap, LinearAlgebra
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

begin
    steel = Material(200e6, 77e6, 7.85, 0.3)
    col = Section(steel, 1.3e-2, 3.0e-4, 1.0e-4, 1.0e-6)
    primary = Section(steel, 1.1e-2, 4.5e-4, 1.2e-4, 1.0e-6)
    joist = Section(steel, 4.0e-3, 6.0e-5, 4.0e-6, 2.0e-7)
    brace = Section(steel, 2.0e-3)
    H = 3.5
    # joistSpacing = 2.0: round(6/2) = 3 joist lines per bay counting both column
    # lines, i.e. ONE joist at midspan of each primary beam (3 m apart)
    bldg = Frame(3, 6.0, 2, 6.0, 3, H, 2.0, col, primary, joist, brace)
    model = bldg.model
    counts = Dict(id => count(el -> el.id == id, model.elements) for id in unique(el.id for el in model.elements))
end

#=
## Load cases

`case = :D / :L / :W` tags each load. Dead: joists 4 kN/m + self-weight of
the beams; live: joists 3 kN/m; wind: 20 kN per storey in +x on the windward
(x = 0) face.
=#
begin
    empty!(model.loads)
    for el in model.elements
        if el.id == :joist
            push!(model.loads, LineLoad(el, [0.0, 0.0, -4.0]; case = :D), LineLoad(el, [0.0, 0.0, -3.0]; case = :L))
        elseif el.id in (:primary, :column)
            push!(model.loads, SelfWeight(el; g = [0.0, 0.0, -9.81], case = :D))
        end
    end
    levels = [k * H for k in 1:3]
    for z in levels
        face = filter(n -> n.position[1] ≈ 0.0 && n.position[3] ≈ z, model.nodes)
        append!(model.loads, [NodeForce(n, [20.0 / length(face), 0.0, 0.0]; case = :W) for n in face])
    end
    cr = solve_cases!(model; cases = [:D, :L, :W])
end

#=
## Storey drift under wind
=#
begin
    rW = cr.results[3]
    drift = [sum(n -> displacement(rW, n)[1], filter(n -> n.position[3] ≈ z, model.nodes)) /
             count(n -> n.position[3] ≈ z, model.nodes) for z in levels]
    interstorey = diff([0.0; drift]) ./ H
    f = fig_bars(["L1", "L2", "L3"], 1000 .* drift; ylabel = "storey displacement [mm]",
        title = "wind drift (max ratio 1/$(round(Int, 1 / maximum(interstorey))))", width = 215)
    savefig("A4_drift.pdf", f)
end

#=
`fig_3d` reads `model.results`; to draw one case, solve a model with that
case's loads (cheap — same nodes and elements):
=#
begin
    windonly = Model(model.nodes, model.elements, filter(l -> l.case == :W, model.loads))
    solve!(windonly)
    f = fig_3d(windonly; scale = 100.0, linecolor = :axial, title = "W: axial force, ×100", width = W_HALF)
    savefig("A4_building_W.pdf", f)
end

#=
## Combinations and an envelope

`LoadCombination`s are exact superpositions of the case results; `envelope`
gives a member's station-wise internal-force extrema over them.
=#
begin
    combos = [LoadCombination(:ULS1, [:D => 1.4]),
              LoadCombination(:ULS2, [:D => 1.2, :L => 1.6]),
              LoadCombination(:ULS3, [:D => 1.2, :L => 1.0, :W => 1.0]),
              LoadCombination(:ULS4, [:D => 0.9, :W => 1.0])]
    girder = first(filter(el -> el.id == :primary, model.elements))
    env = envelope(model, girder, cr, combos; resolution = 40)
    # pick the bending component that the gravity loads activate (rows or columns are force components)
    comp = size(env.lo, 1) == length(env.x) ? axes(env.lo, 2) : axes(env.lo, 1)
    getc(Mat, k) = size(Mat, 1) == length(env.x) ? Mat[:, k] : Mat[k, :]
    k = argmax([maximum(abs, getc(env.hi, k)) + maximum(abs, getc(env.lo, k)) for k in comp if k != 1])
    k = collect(comp)[k + 1]
    f = Figure(size = (W_HALF, round(Int, 0.8W_HALF)))
    ax = Axis(f[1, 1]; aspect = nothing, xlabel = "{x} along girder [m]", ylabel = "moment [kN·m]", title = "moment envelope")
    band!(ax, env.x, getc(env.lo, k), getc(env.hi, k); color = (kjl_black, 0.5))     # one internal force: black, 50 % fill
    lines!(ax, env.x, getc(env.hi, k); color = kjl_black, linewidth = 1); lines!(ax, env.x, getc(env.lo, k); color = kjl_black, linewidth = 1)
    hlines!(ax, [0.0]; color = kjl_midgray, linewidth = 0.5, linestyle = :solid)
    savefig("A4_envelope.pdf", f)
end

#=
## Numbers for the report

`reports/A4_frame_3d_building.pdf` prints these (see `figures.jl`, `report_data`).
Statics check: the support reactions of each case against its total load.
=#
begin
    joists = filter(el -> el.id == :joist, model.elements)
    Ljoist = sum(length, joists)
    sw_total = sum(el -> ρA(el.section) * 9.81 * length(el), filter(el -> el.id in (:primary, :column), model.elements))
    ground = filter(n -> n.position[3] ≈ 0.0, model.nodes)
    Rsum(case, dof) = sum(n -> reaction(cr[case], n)[dof], ground)
    compname = ("N", "Vy", "Mz", "Vz", "My", "T")[k]
    percombo = [(name = string(c.name), factors = [(case = string(p.first), factor = p.second) for p in c.factors],
                 max = maximum(getc(envelope(model, girder, cr, [c]; resolution = 40).hi, k)),
                 min = minimum(getc(envelope(model, girder, cr, [c]; resolution = 40).lo, k))) for c in combos]
    report_data("A4";
        nx = 3, dx = 6.0, ny = 2, dy = 6.0, storeys = 3, H = H, joist_spacing = 2.0,
        plan_x = 3 * 6.0, plan_y = 2 * 6.0, E = steel.E,
        sections = [(name = "column", A = col.A, Ix = col.Ix, Iy = col.Iy),
                    (name = "primary", A = primary.A, Ix = primary.Ix, Iy = primary.Iy),
                    (name = "joist", A = joist.A, Ix = joist.Ix, Iy = joist.Iy),
                    (name = "brace", A = brace.A, Ix = brace.Ix, Iy = brace.Iy)],
        counts = [(id = string(id), n = counts[id]) for id in (:column, :primary, :joist, :brace)],
        n_nodes = length(model.nodes), n_elements = length(model.elements), n_loads = length(model.loads),
        n_supports = length(ground),
        joist_D = 4.0, joist_L = 3.0, wind_per_storey = 20.0,
        joist_length_total = Ljoist,
        totals = (D = 4.0 * Ljoist + sw_total, D_selfweight = sw_total, L = 3.0 * Ljoist, W = 20.0 * length(levels)),
        reactions = (D = Rsum(:D, 3), L = Rsum(:L, 3), W = Rsum(:W, 1)),
        levels = levels, drift_mm = 1000 .* drift, interstorey = interstorey,
        max_interstorey_inv = 1 / maximum(interstorey),
        girder_length = length(girder), component = compname, n_stations = length(env.x),
        envelope_max = maximum(getc(env.hi, k)), envelope_min = minimum(getc(env.lo, k)),
        combos = percombo)
end
