#=
Walkthrough 1 — a shear wall with openings, standalone

Shift+enter through the `begin … end` blocks top to bottom (or run the file).
Units: kN, m.

A two-storey concrete shear wall (6 m long, 2 × 3 m tall, 200 mm thick)
with a ground-floor door and a first-floor window. It carries its own
weight, floor gravity at each level, and a lateral storey force at each
level. We mesh it, solve it, and read it the way a designer does: drift,
membrane forces, compression struts, and pier forces at a section cut.

Environment: `examples/Project.toml` (Asap and AsapOptim by path + plotting).
    julia --project=examples examples/shells/01_wall_basics.jl
=#

begin
    using Asap, LinearAlgebra
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

#=
## Material and section

`ShellSection(material, t)` is a homogeneous isotropic shell. For cracked
concrete use `ShellRigiditySection(section; m11 = …, f11 = …)` — the
property-modifier convention of building-analysis software (which values a
code prescribes is your call, not Asap's).
=#
begin
    concrete = Material(30.0e6, 2.4, 0.2)                 # E [kN/m²], ρ [t/m³], ν
    wall = ShellSection(concrete, 0.20)
    wall_cracked = ShellRigiditySection(wall; f11 = 0.7, f22 = 0.7, f12 = 0.7)   # e.g. a cracked in-plane stiffness study
    wall
end

#=
## Mesh

`mesh_wall(origin, along, up, length, height, section; …)` builds a
structured mesh. Openings are rectangles `(s₀, z₀, width, height)` in wall
coordinates; their edges, the storey level `z_lines`, and any `s_lines`
become grid lines, so piers and spandrels are meshed exactly.
=#
begin
    L, H1, H2 = 6.0, 3.0, 3.0
    door = (1.0, 0.0, 1.0, 2.2)
    window = (3.0, 4.0, 1.5, 1.2)
    w = mesh_wall([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], L, H1 + H2, wall;
        openings = [door, window], z_lines = [1.0, H1], size = 0.25)
end

#=
## Loads

- self-weight of every shell (`SelfWeight`, from ρ·t),
- floor gravity: 30 kN/m along each level, lumped to the level's nodes,
- lateral storey forces: 120 kN at level 1, 160 kN at the roof, applied
  uniformly along each level (a rigid diaphragm would do this for a whole
  building — see walkthrough 3).
=#
begin
    level(z) = [w.grid[i, findfirst(==(z), w.z)] for i in eachindex(w.s) if w.grid[i, findfirst(==(z), w.z)] !== nothing]
    lvl1, roof = level(H1), level(H1 + H2)
    loads = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0.0, -9.81]) for el in w.elements]...]
    for (nodes, V) in ((lvl1, 120.0), (roof, 160.0))
        n = length(nodes)
        append!(loads, [NodeForce(nd, [V / n, 0.0, -30.0 * L / n]) for nd in nodes])
    end
    model = Model(w.nodes, AbstractElement{Float64}[w.elements...], loads)
    solve!(model)
    mean_x(nodes) = sum(n -> displacement(model.results, n)[1], nodes) / length(nodes)
    drifts = (level1_drift = mean_x(lvl1), roof_drift = mean_x(roof), roof_drift_ratio = mean_x(roof) / (H1 + H2))
    @show drifts
end

#=
## Membrane forces and the deformed shape

`ShellGeo(model)` extracts plot-ready arrays: centre resultants in each
element's local axes (e₁ along the wall, e₂ up the wall here), principal
forces, displacements. Nyy is the vertical membrane force: compression at
the leeward toe, tension at the windward heel, concentrations at opening
corners.
=#
begin
    g = ShellGeo(model)
    f = fig_wall(model, w; geo = g, field = :Nyy, scale = 200.0,
        title = "vertical membrane force Nyy, deformed ×200")
    savefig("01_wall_Nyy.png", f)
end

#=
## Compression struts

Principal membrane-force trajectories are the strut-and-tie reading of the
wall: blue segments (compression) are struts, pink (tension) the ties that
reinforcement must carry — note the struts bending around the openings and
the tension field below the window.
=#
begin
    f = fig_principal(model, w; geo = g)
    savefig("01_wall_principal.png", f)
end

#=
## Pier forces by section cut

`wall_cut(w, z; s_range)` is the free body below height `z` within
`s_range` (a pier between openings), and `section_cut` integrates the
forces crossing the cut from element nodal forces — exact statics on any
mesh, no stress smoothing. Here at z = 1 m, through the door: the left pier
(0–1 m) and the right pier (2–6 m).
=#
begin
    eL, cL = wall_cut(w, 1.0; s_range = (0.0, 1.0))
    eR, cR = wall_cut(w, 1.0; s_range = (2.0, L))
    # report in pier axes: rows = (along wall, through thickness, vertical)
    axes_pier = [1.0 0 0; 0 1 0; 0 0 1]
    FL, ML = section_cut(model, eL, cL; axes = axes_pier)
    FR, MR = section_cut(model, eR, cR; axes = axes_pier)
    piers = (left = (V = FL[1], P = FL[3], M = ML[2]), right = (V = FR[1], P = FR[3], M = MR[2]))
    @show FL[1] + FR[1]          # = total storey shear above z = 1 (280 kN), exactly
    f = fig_bars(["left pier V", "right pier V"], [FL[1], FR[1]]; ylabel = "shear [kN]",
        title = "pier shears at z = 1 m (sum = 280 kN)")
    savefig("01_wall_pier_shears.png", f)
    @show piers
end

#=
## Resultants at any point of an element

`shell_resultants(model, el; at = (ξ, η))` — natural coordinates in
[−1, 1]²; centres are the superconvergent sampling points.
=#
begin
    corner_el = only(filter(el -> any(n -> n.position ≈ [2.0, 0.0, 2.2], el.nodes) &&
                                 all(n -> n.position[1] >= 2.0 && n.position[3] <= 2.2, el.nodes), w.elements))
    shell_resultants(model, corner_el; at = (-1.0, 1.0))     # at the door's top-right corner
end

#=
## Mesh convergence

The drilling membrane is accurate on coarse meshes — a wall meshed two
elements wide is already close. Roof drift against element size:
=#
begin
    sizes = [1.0, 0.5, 0.25, 0.125]
    drift = map(sizes) do h
        wh = mesh_wall([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], L, H1 + H2, wall;
            openings = [door, window], z_lines = [1.0, H1], size = h)
        lv(z) = [wh.grid[i, findfirst(==(z), wh.z)] for i in eachindex(wh.s) if wh.grid[i, findfirst(==(z), wh.z)] !== nothing]
        ld = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0.0, -9.81]) for el in wh.elements]...]
        for (nodes, V) in ((lv(H1), 120.0), (lv(H1 + H2), 160.0))
            append!(ld, [NodeForce(nd, [V / length(nodes), 0.0, -30.0 * L / length(nodes)]) for nd in nodes])
        end
        mh = Model(wh.nodes, AbstractElement{Float64}[wh.elements...], ld)
        solve!(mh)
        sum(n -> displacement(mh.results, n)[1], lv(H1 + H2)) / length(lv(H1 + H2))
    end
    f = fig_series(sizes, [("roof drift", 1000 .* drift)]; xlabel = "element size [m]",
        ylabel = "roof drift [mm]", xscale = log2, title = "mesh convergence")
    savefig("01_wall_convergence.png", f)
    drift
end

#=
## Numbers for the report

`reports/S1_wall_basics.pdf` prints these (see `figures.jl`, `report_data`).
=#
begin
    corner = shell_resultants(model, corner_el; at = (-1.0, 1.0))
    sw = areal_mass(wall) * 9.81                     # self-weight per wall area [kN/m²]
    net_area = L * (H1 + H2) - door[3] * door[4] - window[3] * window[4]
    above_cut = L * (H1 + H2 - 1.0) - door[3] * (door[4] - 1.0) - window[3] * window[4]
    report_data("S1";
        length = L, H1 = H1, H2 = H2, thickness = wall.t,
        E = concrete.E, nu = concrete.ν, rho = concrete.ρ, cracked_factor = 0.7,
        door = door, window = window, mesh_size = 0.25,
        n_elements = length(w.elements), n_nodes = length(w.nodes),
        n_level1 = length(lvl1), n_roof = length(roof),
        floor_load = 30.0, storey_forces = (level1 = 120.0, roof = 160.0),
        self_weight_per_area = sw, self_weight_total = sw * net_area,
        drifts_mm = (level1 = 1000 * drifts.level1_drift, roof = 1000 * drifts.roof_drift),
        roof_drift_ratio = drifts.roof_drift_ratio,
        cut_z = 1.0,
        piers = (left = piers.left, right = piers.right),
        pier_shear_sum = FL[1] + FR[1], storey_shear_above_cut = 120.0 + 160.0,
        pier_axial_sum = FL[3] + FR[3], gravity_above_cut = sw * above_cut + 2 * 30.0 * L,
        corner = (N = corner.N, M = corner.M, Q = corner.Q),
        convergence = (sizes = sizes, roof_drift_mm = 1000 .* drift))
end
