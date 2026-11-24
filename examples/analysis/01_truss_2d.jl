#=
Walkthrough A1 — a planar truss, by hand and by generator

Shift+enter through the `begin … end` blocks. Units: kN, m.

A 12 m Pratt bridge truss, 2 m deep: nodes, supports, axial-only members,
deck loads, axial forces and reactions, the deflected shape — then the same
idea from a one-line generator, and a quick depth study.
=#

begin
    using Asap, LinearAlgebra
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

#=
## Material and section

`Section(material, A)` is an axial-only section — all a `TrussElement` uses.
=#
begin
    steel = Material(200e6, 77e6, 7.85, 0.3)        # E [kN/m²], G [kN/m²], ρ [t/m³], ν
    bar = Section(steel, 2.0e-3)                     # 20 cm²
end

#=
## Geometry, supports, members

Nodes carry their supports as 6 Booleans (true = free): the left support is
pinned, the right a roller. `planarize!` then fixes every node's out-of-plane
DOFs, so the truss behaves as 2D in the XY plane. Truss nodes' rotations
are retired automatically (nothing stiffens them).
=#
begin
    L, H, n = 12.0, 2.0, 6
    bottom = [Node([L * i / n, 0.0, 0.0], :free, :bottom) for i in 0:n]
    top = [Node([L * i / n, H, 0.0], :free, :top) for i in 1:n-1]
    bottom[1] = Node([0.0, 0.0, 0.0], Bool[0, 0, 0, 1, 1, 1], :pin)        # pinned
    bottom[end] = Node([L, 0.0, 0.0], Bool[1, 0, 0, 1, 1, 1], :roller)     # vertical roller
    nodes = vcat(bottom, top)
    planarize!(nodes, :XY)

    T(a, b, id) = TrussElement(a, b, bar, id)
    members = AbstractElement{Float64}[]
    append!(members, [T(bottom[i], bottom[i+1], :bottomchord) for i in 1:n])
    append!(members, [T(top[i], top[i+1], :topchord) for i in 1:n-2])
    push!(members, T(bottom[1], top[1], :endpost), T(top[end], bottom[end], :endpost))
    append!(members, [T(bottom[i+1], top[i], :vertical) for i in 1:n-1])
    # Pratt diagonals: sloping down towards midspan (tension under gravity)
    append!(members, [T(top[i], bottom[i+2], :diagonal) for i in 1:n÷2-1])
    append!(members, [T(top[i+1], bottom[i+1], :diagonal) for i in n÷2:n-2])
    length(members)
end

#=
## Loads and solve

20 kN at every interior bottom-chord node (the deck).
=#
begin
    loads = AbstractLoad{Float64}[NodeForce(nd, [0.0, -20.0, 0.0]) for nd in bottom[2:end-1]]
    truss = Model(nodes, members, loads)
    solve!(truss)
    res = truss.results
    (left_reaction = reaction(res, bottom[1])[2], right_reaction = reaction(res, bottom[end])[2],
     midspan_deflection_mm = 1000 * displacement(res, bottom[n÷2+1])[2])
end

#=
## Axial forces

`axial_force(results, element)` — tension positive. Statics check (method of
sections): the bottom chord in the panel x = 4–6 m, cut together with the
diagonal and top chord meeting at the top node x = 4 m, carries
M(4 m)/H = (50·4 − 20·2)/2 = 80 kN.
=#
begin
    midchord = members[n÷2]
    @show axial_force(res, midchord)
    byid = Dict(id => extrema(axial_force(res, el) for el in members if el.id == id)
                for id in unique(el.id for el in members))
end

begin
    f = fig_lines(truss; plane = :XY, scale = 200.0, labels = true,
        title = "Pratt truss: axial force [kN], deflection ×200")
    savefig("A1_pratt.png", f)
end

#=
## The same idea from a generator

`Warren2D(n, dx, dy, section; load)` builds, loads and solves a Warren
truss in one line (see the other generators: `Pratt2D`, `SpaceFrame`,
`Frame`, ground structures, …).
=#
begin
    warren = Warren2D(11, 1.0, 1.5, bar; load = [0.0, -15.0, 0.0])
    f = fig_lines(warren.model; plane = :XY, scale = 100.0, title = "Warren2D generator, deflection ×100")
    savefig("A1_warren.png", f)
end

#=
## Depth study

Deeper trusses have smaller chord forces but longer webs — midspan
deflection against depth, same members and loads.
=#
begin
    function pratt_deflection(H)
        b = [Node([L * i / n, 0.0, 0.0], :free) for i in 0:n]
        t = [Node([L * i / n, H, 0.0], :free) for i in 1:n-1]
        b[1] = Node([0.0, 0.0, 0.0], Bool[0, 0, 0, 1, 1, 1]); b[end] = Node([L, 0.0, 0.0], Bool[1, 0, 0, 1, 1, 1])
        nd = vcat(b, t); planarize!(nd, :XY)
        els = AbstractElement{Float64}[[TrussElement(b[i], b[i+1], bar) for i in 1:n]...,
            [TrussElement(t[i], t[i+1], bar) for i in 1:n-2]..., TrussElement(b[1], t[1], bar), TrussElement(t[end], b[end], bar),
            [TrussElement(b[i+1], t[i], bar) for i in 1:n-1]...,
            [TrussElement(t[i], b[i+2], bar) for i in 1:n÷2-1]..., [TrussElement(t[i+1], b[i+1], bar) for i in n÷2:n-2]...]
        m = Model(nd, els, AbstractLoad{Float64}[NodeForce(x, [0.0, -20.0, 0.0]) for x in b[2:end-1]])
        solve!(m)
        return -1000 * displacement(m.results, b[n÷2+1])[2]
    end
    depths = 1.0:0.25:3.0
    f = fig_series(collect(depths), [("midspan deflection", pratt_deflection.(depths))];
        xlabel = "depth H [m]", ylabel = "deflection [mm]", title = "Pratt truss: depth study", legendpos = :rt)
    savefig("A1_depth.png", f)
end

#=
## Numbers for the report

`reports/A1_truss_2d.typ` prints these (see `figures.jl`, `report_data`).
=#
begin
    report_data("A1";
        span = L, depth = H, panels = n, n_members = length(members), panel_load = 20.0,
        E = steel.E, area = bar.A,
        reactions = (left = reaction(res, bottom[1])[2], right = reaction(res, bottom[end])[2]),
        midspan_deflection_mm = -1000 * displacement(res, bottom[n÷2+1])[2],
        midchord_force = axial_force(res, midchord), midchord_hand = (50 * 4 - 20 * 2) / 2,
        forces = [(group = string(id), min = byid[id][1], max = byid[id][2]) for id in (:bottomchord, :topchord, :endpost, :vertical, :diagonal)],
        warren = (bays = 11, n_members = length(warren.model.elements),
                  max_deflection_mm = 1000 * maximum(nd -> abs(displacement(warren.model.results, nd)[2]), warren.model.nodes)),
        depth_study = (depths = collect(depths), deflection_mm = pratt_deflection.(depths)))
end
