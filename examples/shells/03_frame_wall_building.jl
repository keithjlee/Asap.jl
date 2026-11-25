#=
Walkthrough 3 — shells WITH frames: a small building with a coupled shear wall

Shift+enter through the blocks. Units: kN, m.

A two-storey building, 12 m (x) × 8 m (y), storeys 3.5 m:

- perimeter FRAMES (columns + beams) on all four sides,
- a COUPLED SHEAR WALL on the centreline x = 6 m: two 2.5 m wall piers with
  a 1 m door between them, joined at each floor by a frame COUPLING BEAM
  that lies in the wall's plane — it connects to the wall at a single node
  at each end, so its end moment can only pass into the wall through the
  shell's DRILLING rotation,
- a RIGID DIAPHRAGM at each floor tying everything together in plan.

Two lateral cases (X along the long frames, Y along the wall), gravity, and
the questions a designer asks: drift, how the base shear splits between the
wall and the frames, pier forces, coupling-beam forces.
=#

begin
    using Asap, LinearAlgebra
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

#=
## Materials and sections

For a readable walkthrough the frame members use equal bending stiffness
about both axes (no roll-angle bookkeeping); real members would set the
section orientation.
=#
begin
    conc = Material(30.0e6, 2.4, 0.2)
    col = Section(conc, 0.16, 2.13e-3, 2.13e-3, 3.6e-3)      # 400 × 400 column
    beam = Section(conc, 0.15, 2.0e-3, 2.0e-3, 2.8e-3)       # 300 × 500 beam (equalized)
    cbeam = Section(conc, 0.12, 1.6e-3, 1.6e-3, 1.8e-3)      # coupling beam
    wallsec = ShellSection(conc, 0.20)
    H = 3.5
    levels = [H, 2H]
end

#=
## The coupled wall

Two `mesh_wall`s in the plane x = 6 (`along` = +y, `up` = +z), each 2.5 m
long, with the storey levels as hard lines so floor nodes exist there.
=#
begin
    pierA = mesh_wall([6.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 2.5, 2H, wallsec;
        z_lines = levels, size = 0.5, id = :pierA)
    pierB = mesh_wall([6.0, 4.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 2.5, 2H, wallsec;
        z_lines = levels, size = 0.5, id = :pierB)
    wallnodes = vcat(pierA.nodes, pierB.nodes)
    wallels = vcat(pierA.elements, pierB.elements)
    # coupling beams: in the wall plane, from pier A's inner edge to pier B's, at each level
    atlevel(w, s, z) = only(filter(n -> n.position ≈ w.origin + s * w.along + z * w.up, w.nodes))
    cbeams = [FrameElement(atlevel(pierA, 2.5, z), atlevel(pierB, 0.0, z), cbeam, :coupling) for z in levels]
    length(wallels), length(cbeams)
end

#=
## The frames

Columns at the perimeter grid, beams between them at each level.
=#
begin
    colxy = [(0.0, 0.0), (4.0, 0.0), (8.0, 0.0), (12.0, 0.0), (0.0, 8.0), (4.0, 8.0), (8.0, 8.0), (12.0, 8.0),
             (0.0, 4.0), (12.0, 4.0)]
    nodeat = Dict{Tuple{Float64,Float64,Float64},Node{Float64}}()
    getnode(x, y, z) = get!(nodeat, (x, y, z)) do
        Node([x, y, z], z == 0 ? :fixed : :free)
    end
    frames = FrameElement[]
    for (x, y) in colxy, (z0, z1) in ((0.0, H), (H, 2H))
        push!(frames, FrameElement(getnode(x, y, z0), getnode(x, y, z1), col, :column))
    end
    beamlines = [[(0.0, 0.0), (4.0, 0.0), (8.0, 0.0), (12.0, 0.0)], [(0.0, 8.0), (4.0, 8.0), (8.0, 8.0), (12.0, 8.0)],
                 [(0.0, 0.0), (0.0, 4.0), (0.0, 8.0)], [(12.0, 0.0), (12.0, 4.0), (12.0, 8.0)]]
    for z in levels, line in beamlines, k in 1:length(line)-1
        (x0, y0), (x1, y1) = line[k], line[k+1]
        push!(frames, FrameElement(getnode(x0, y0, z), getnode(x1, y1, z), beam, :beam))
    end
    framenodes = collect(values(nodeat))
    length(frames)
end

#=
## Rigid diaphragms

One per floor: slaves are every node at that level (column tops and the
wall's floor-line nodes). `RigidDiaphragm` creates its master at the
centroid; storey forces and accidental torsion go there.
=#
begin
    allnodes = vcat(wallnodes, framenodes)
    floors = [RigidDiaphragm(filter(n -> n.position[3] ≈ z, allnodes); id = Symbol("L", i)) for (i, z) in enumerate(levels)]
end

#=
## Loads (as cases)

- `:D` — self-weight of walls and frames, and floor dead load as a line
  load on the perimeter beams (15 kN/m);
- `:EX`, `:EY` — storey forces 150 kN (level 1) and 250 kN (roof) at the
  diaphragm masters, plus 5% accidental-eccentricity torsion.
=#
begin
    g = [0.0, 0.0, -9.81]
    elements = AbstractElement{Float64}[wallels..., cbeams..., frames...]
    loads = AbstractLoad{Float64}[
        [SelfWeight(el; g = g, case = :D) for el in wallels]...,
        [SelfWeight(el; g = g, case = :D) for el in frames]...,
        [LineLoad(el, [0.0, 0.0, -15.0]; case = :D) for el in frames if el.id == :beam]...]
    for (d, V) in zip(floors, (150.0, 250.0))
        push!(loads, NodeForce(d.master, [V, 0, 0]; case = :EX), NodeMoment(d.master, [0, 0, 0.05 * 8.0 * V]; case = :EX))
        push!(loads, NodeForce(d.master, [0, V, 0]; case = :EY), NodeMoment(d.master, [0, 0, 0.05 * 12.0 * V]; case = :EY))
    end
    model = Model(allnodes, elements, loads; diaphragms = floors)
    cr = solve_cases!(model; cases = [:D, :EX, :EY])
    model
end

#=
## Drift

Average storey displacement = the master's translation. The Y direction
(wall) is far stiffer than X (frames only).
=#
begin
    rD, rEX, rEY = cr.results
    drift(res, dir) = [displacement(res, d.master)[dir] for d in floors]
    drifts = (EX_mm = 1000 .* drift(rEX, 1), EY_mm = 1000 .* drift(rEY, 2),
              EX_ratio_roof = drift(rEX, 1)[2] / 2H, EY_ratio_roof = drift(rEY, 2)[2] / 2H)
    @show drifts
end

#=
## Who carries the Y shear — the wall or the frames?

Base reactions split by what they support: the wall piers' base nodes vs
the columns' bases. (Asap's reaction convention: elastic forces at the
supports.)
=#
begin
    basewall = filter(n -> n.position[3] == 0, wallnodes)
    basecols = filter(n -> n.position[3] == 0, framenodes)
    Vwall = -sum(n -> reaction(rEY, n)[2], basewall)
    Vframe = -sum(n -> reaction(rEY, n)[2], basecols)
    share = (wall = Vwall, frames = Vframe, total = Vwall + Vframe, wall_fraction = Vwall / (Vwall + Vframe))
    @show share
    f = fig_bars(["coupled wall", "frames"], [Vwall, Vframe]; ylabel = "base shear [kN]",
        title = "EY: $(round(100 * share.wall_fraction, digits = 1))% to wall", highlight = 1, width = 180)
    savefig("03_building_shear_share.pdf", f)
end

#=
## Pier forces (EY) at the first floor

`wall_cut` + `section_cut` on each pier, reported in pier axes (along the
wall, through it, vertical): the two piers work as a coupled pair — one in
tension, one in compression, with the couple carried through the coupling
beams.
=#
begin
    pieraxes = [0.0 1 0; 1 0 0; 0 0 1]           # rows: along wall (y), through (x), vertical (z)
    function pierforce(w, res)
        els, cut = wall_cut(w, H)
        F, M = section_cut(model, els, cut; axes = pieraxes, results = res)
        (V = F[1], P = F[3], M = M[2])
    end
    piers = (A = pierforce(pierA, rEY), B = pierforce(pierB, rEY))
    @show piers
end

#=
## Coupling beams

A coupling beam is a frame element connected to each pier at ONE node.
Its end moments (local My/Mz) are real — they pass into the wall through
the shell's drilling rotation. A membrane without a drilling DOF would
leave the beam pinned to the wall and the piers uncoupled.
=#
begin
    cb = [element_forces(rEY, b) for b in cbeams]
    coupling = [(level = i, shear = maximum(abs, f[[2, 3]]), end_moment = maximum(abs, f[[5, 6, 11, 12]])) for (i, f) in enumerate(cb)]
    @show coupling
end

#=
## Figures — the building under EY, and the wall's force flow
=#
begin
    rg = ShellGeo(model; results = rEY)
    f = fig_3d(model; geo = rg, field = :umag, scale = 300.0, title = "EY, deformed ×300", width = 312)
    savefig("03_building_EY_3d.pdf", f)
end

begin
    # the full wall line (both piers) in one elevation: build a geo once, plot each pier
    f = fig_wall(model, pierA; geo = rg, field = :Nyy, scale = 300.0, title = "pier A: N_{yy}, EY", width = 196)
    savefig("03_pierA_Nyy_EY.pdf", f)
end

begin
    f = fig_principal(model, pierB; geo = rg, title = "pier B: principal, EY", width = 196)
    savefig("03_pierB_principal_EY.pdf", f)
end

#=
## Combinations

Results superpose: `combine` gives any factored combination exactly, and
`shell_envelope` the resultant extrema of an element over combinations.
=#
begin
    combos = [LoadCombination(:ULS1, [:D => 1.2, :EY => 1.0]), LoadCombination(:ULS2, [:D => 0.9, :EY => -1.0])]
    heel = pierA.elements[1]                       # bottom-left element of pier A
    env = shell_envelope(model, heel, cr, combos)
    @show env.Nmin env.Nmax
end

#=
## Numbers for the report

`reports/S3_frame_wall_building.pdf` prints these (see `figures.jl`, `report_data`).
=#
begin
    secinfo(s) = (A = s.A, Ix = s.Ix, Iy = s.Iy, J = s.J)
    report_data("S3";
        plan_x = 12.0, plan_y = 8.0, wall_x = 6.0, column_spacing = 4.0, storey_height = H, n_storeys = length(levels),
        E = conc.E, nu = conc.ν, rho = conc.ρ,
        column = secinfo(col), beam = secinfo(beam), coupling_beam = secinfo(cbeam),
        wall_thickness = wallsec.t, pier_length = 2.5, door_width = 1.0, wall_mesh_size = 0.5,
        n_wall_elements = length(wallels), n_frames = length(frames), n_columns = 2 * length(colxy),
        n_coupling = length(cbeams), n_slaves = [length(d.nodes) for d in floors],
        beam_line_load = 15.0, storey_forces = (level1 = 150.0, roof = 250.0), eccentricity = 0.05,
        torsion_EX = [0.05 * 8.0 * V for V in (150.0, 250.0)], torsion_EY = [0.05 * 12.0 * V for V in (150.0, 250.0)],
        drifts = (EX_mm = drifts.EX_mm, EY_mm = drifts.EY_mm,
                  EX_ratio_roof = drifts.EX_ratio_roof, EY_ratio_roof = drifts.EY_ratio_roof),
        stiffness_ratio = drifts.EX_mm[2] / drifts.EY_mm[2],
        share = share, applied_EY = 150.0 + 250.0,
        piers = (A = piers.A, B = piers.B),
        coupling = [(level = c.level, shear = c.shear, end_moment = c.end_moment) for c in coupling],
        coupling_shear_sum = sum(c.shear for c in coupling),
        envelope = (combos = ["1.2 D + 1.0 EY", "0.9 D − 1.0 EY"],
                    Nmin = env.Nmin, Nmax = env.Nmax, Mmin = env.Mmin, Mmax = env.Mmax))
end
