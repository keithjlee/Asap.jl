# Shells, walls and diaphragms

```@setup shells
using Asap
```

Asap's [`ShellElement`](@ref) is a flat shell in two shapes, which mix freely in one mesh:

- **Quadrilateral** (4 nodes): a drilling-rotation membrane (Ibrahimbegović–Taylor–Wilson, with two complementary incompatible modes) combined with the MITC4 Reissner–Mindlin plate.
- **Triangle** (3 nodes): Felippa's ANDES-OPT drilling membrane combined with the MITC3+ Reissner–Mindlin plate. ANDES-OPT is optimal in energy for in-plane bending of rectangular two-triangle units (the strain energy is exact). The pointwise stresses are not exact: the centroid membrane force of a two-triangle-per-cell beam is off by about 1/(6·rows) of the extreme-fibre value.

In a mixed mesh, a triangle uses the quadrilateral's drilling force lumping on every edge it shares with a quadrilateral (the element's `drilling` field, set by [`process!`](@ref)). Both membranes then lump a constant stress state alike, so mixed meshes pass the patch test exactly. Only those transition triangles give up the ANDES bending optimality.

Shells mix freely with frames and trusses, and the drilling DOF lets an in-plane beam transfer moment into a wall at a shared node.

The elements are verified on patch tests, the MacNeal–Harder benchmarks, Timoshenko & Woinowsky-Krieger plates and the Scordelis–Lo roof; the test suite runs every one of them.

## Sections

A [`ShellSection`](@ref) is a homogeneous isotropic shell of thickness `t`. A [`ShellRigiditySection`](@ref) stores the thickness-integrated rigidities directly. Use it for cracked or effective concrete and for orthotropic equivalent plates. Its modifier constructor applies the property-modifier convention of building-analysis software:

```@example shells
concrete = Material(3.0e7, 2.4, 0.2)           # kN, m
wall = ShellSection(concrete, 0.2)
cracked = ShellRigiditySection(wall; m11 = 0.35, m22 = 0.35, m12 = 0.35)
```

## Meshing a wall

[`mesh_wall`](@ref) builds a structured mesh of a planar rectangular wall:
- rectangular openings are removed;
- storey levels, column lines and opening edges are guaranteed grid lines;
- nodes passed in `reuse` are shared, which is how walls join each other and frames.

```@example shells
w = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 6.0, 3.0, wall;
              openings = [(2.0, 0.0, 1.0, 2.1)],      # a door: s₀, z₀, width, height
              z_lines = [1.0], size = 0.5)
```

## General planar meshes

For any other outline, [`mesh_polygon`](@ref) triangulates a polygon with polygonal holes (constrained Delaunay with quality refinement; load DelaunayTriangulation.jl to enable it). Put a vertex wherever a support or a load starts or ends (or pass the point in `points`: a point on an edge splits that edge, so it stays on it even on a sloped edge, where rounding puts it a hair off the line), and use [`circle_polygon`](@ref) for round openings. Meshes from any other mesher (points and triangle/quad cells) become elements through [`shell_mesh`](@ref), and Gmsh files (MSH 4.1 ASCII) through [`read_gmsh`](@ref), which keeps Gmsh's physical groups as named node and element sets ([`group_nodes`](@ref), [`group_elements`](@ref)) and reads curved surfaces as spatial meshes. `mesh_polygon` is reproducible: the same input and `seed` give the same mesh. The target edge length `h` is a number or a function `(x, y) -> h` of the in-plane position, which grades the mesh: fine near loads, bearings and re-entrant corners, coarse elsewhere. Both return a [`ShellMesh`](@ref), which keeps the in-plane coordinates for region queries:

- [`nodes_on_segment`](@ref) and [`boundary_nodes`](@ref): nodes along a bearing or a loaded edge;
- [`edge_loads`](@ref) with a segment: consistent line loads along a boundary, drilling moments included;
- [`elements_in`](@ref): elements in a region (a bearing plate's patch load, a zone to read stresses from);
- [`total_reaction`](@ref): the resultant of the reactions over a node group, e.g. one bearing.

To integrate the force across a line through the material (a nib/beam interface, the level a hanger tie crosses), pass the line to `mesh_polygon` as one of its `segments`: element edges then follow it, and [`section_cut`](@ref) with the elements on one side and [`nodes_on_segment`](@ref) gives the exact force and moment across it. Summing the cut one node at a time along the line gives the force profile across it, whose tension part is a strut-and-tie model's tie (walkthrough 7, a dapped end).

A flexible bearing (an elastomeric pad) presses with a nearly uniform pressure and cannot pull. Holding the pad's nodes fixed makes it rigid and bilateral instead. To model the pad, apply the reaction as an `edge_loads` pressure and hold the panel only against rigid-body motion (one pin, one roller): the load set balances by itself, and the supports' reactions measure how exactly.

```julia
using DelaunayTriangulation      # activates mesh_polygon
outline = [(0.0, 0.0), (0.5, 0.0), (5.5, 0.0), (6.0, 0.0), (6.0, 3.0), (0.0, 3.0)]
panel = mesh_polygon(outline, wall; holes = [circle_polygon((3.0, 1.5), 0.4)], h = 0.25,
                     e1 = [1.0, 0, 0], e2 = [0.0, 0, 1])          # a panel in the XZ plane
for n in panel.nodes
    fixnode!(n, Bool[1, 0, 1, 0, 1, 0])                        # held out of plane
end
# bearings: translations held along the segment, so the drilling rotation (about Y) too
foreach(n -> fixnode!(n, Bool[0, 0, 0, 0, 0, 0]), nodes_on_segment(panel, (0.0, 0.0), (0.5, 0.0)))
foreach(n -> fixnode!(n, Bool[1, 0, 0, 0, 0, 0]), nodes_on_segment(panel, (5.5, 0.0), (6.0, 0.0)))
loads = AbstractLoad{Float64}[edge_loads(panel, (0.0, 3.0), (6.0, 3.0), [0.0, 0.0, -40.0])...]
```

Where translations are held along a segment of several nodes (a bearing, a clamped edge), hold the drilling rotation there too: the material cannot rotate along a held line, and leaving the drilling DOF free leaves the edge's drilling moments unbalanced. A single-node pin keeps it free.

## Loads

- [`AreaLoad`](@ref) applies a uniform surface traction, in global axes or as a pressure along the element normal.
- `SelfWeight` applies to shells too (`ρt·g`).
- Both lower to consistent nodal loads, including the work in-plane traction does on the drilling rotations.
- For line tractions along an element edge, use [`edge_loads`](@ref).

```@example shells
top = [w.grid[i, end] for i in eachindex(w.s)]
loads = AbstractLoad{Float64}[
    [SelfWeight(el; g = [0.0, 0, -9.81]) for el in w.elements]...,
    [NodeForce(n, [60.0 / length(top), 0, 0]) for n in top]...]
model = Model(w.nodes, AbstractElement{Float64}[w.elements...], loads)
solve!(model)
maximum(abs, model.results.u)
```

## Results: resultants and section cuts

[`shell_resultants`](@ref) returns the membrane forces `N`, moments `M` and transverse shears `Q` (per unit width, local axes) at any point of an element.

[`section_cut`](@ref) integrates the forces transmitted across a cut from element nodal forces, which makes it exact for statics on any mesh. [`wall_cut`](@ref) builds the free body for a pier between openings:

```@example shells
eL, cL = wall_cut(w, 1.0; s_range = (0.0, 2.0))      # left pier at z = 1
eR, cR = wall_cut(w, 1.0; s_range = (3.0, 6.0))      # right pier
FL, ML = section_cut(model, eL, cL)
FR, MR = section_cut(model, eR, cR)
(left_shear = FL[1], right_shear = FR[1], total = FL[1] + FR[1])
```

```@example shells
shell_resultants(model, w.elements[1]).N               # (Nxx, Nyy, Nxy)
```

## Precast joints: links

A [`LinkElement`](@ref) is a zero-length 6-spring connection between coincident nodes, for horizontal or vertical panel joints, bearing pads and dowels. Mesh each panel separately, with `fixity_bottom = :free` on the upper panel so its joint nodes are new, and connect the duplicate nodes along the joint with links.

## Rigid diaphragms

A [`RigidDiaphragm`](@ref) ties the in-plane translations and the rotation about the plane normal of its slave nodes to a master node, imposed exactly by transformation. Apply storey forces and torques to the master:

```@example shells
mat = Material(2.0e8, 7.85, 0.3)
base = [Node([x, y, 0.0], :fixed) for (x, y) in ((0.0, 0.0), (6.0, 0.0), (0.0, 4.0), (6.0, 4.0))]
tops = [Node([x, y, 3.0], :free) for (x, y) in ((0.0, 0.0), (6.0, 0.0), (0.0, 4.0), (6.0, 4.0))]
cols = AbstractElement{Float64}[FrameElement(base[k], tops[k], Section(mat, 1e-2, 1e-4k, 1e-4k, 5e-5)) for k in 1:4]
floor = RigidDiaphragm(tops)
storey = Model(vcat(base, tops), cols,
               AbstractLoad{Float64}[NodeForce(floor.master, [10.0, 0, 0])]; diaphragms = [floor])
solve!(storey)
displacement(storey.results, floor.master)[[1, 2, 6]]   # ux, uy, θz: torsion from the stiffness eccentricity
```

## Envelopes and plot extraction

- [`shell_envelope`](@ref) gives the component-wise resultant extrema of an element over load
  combinations (the shell counterpart of [`envelope`](@ref)).
- [`ShellGeo`](@ref) extracts plot-ready arrays: quads, displacements, centre resultants and
  principal forces ([`principal_forces`](@ref)).
- [`nodal_average`](@ref) smooths any per-element field to the nodes for contour plots, with equal, area or angle weights (`weight = :area` on graded meshes). Average the Cartesian components and derive principal values from them, never the principal values or angles themselves.
- [`Geo`](@ref) / `ModelGeo` cover the frame and truss elements of a mixed model.

## Walkthroughs

`examples/shells/` in the repository has shift+enter walkthroughs with figure generators:
1. a shear wall with openings;
2. slabs;
3. a frame + coupled-wall building with rigid diaphragms;
4. precast joints;
5. thickness sensitivities and optimization.

## Differentiable analysis

Shells take part in the pure path: [`compliance`](@ref)/[`solve`](@ref) on a [`ModelState`](@ref) differentiate with respect to node positions and shell sections (thickness), with shell surface loads re-evaluated per state (self-weight follows `ρt`, pressure follows the element area and orientation).

For fast reverse mode, load ForwardDiff alongside Zygote. The shell kernels' rules (`AsapShellADExt`) evaluate each element's derivative with ForwardDiff instead of tracing it, many times faster. Mooncake (`AsapShellMooncakeExt`) and Enzyme reverse mode (`AsapShellEnzymeExt`) use the same rules; both fail on the shell kernels natively. All backends agree to 1e-16; Enzyme is the fastest, with Zygote close behind (see `examples/ad/ad_backends_shells.jl`).

In AsapOptim, `ThicknessVariable(elements, value, lb, ub)` makes a wall or slab thickness a design variable, and `shell_volume(x, p)` gives the material volume. `solution_tangents` (implicit differentiation) supports shells and rigid diaphragms.

```julia
using Zygote, ForwardDiff
st = extract_state(model)
idx = [findfirst(e -> e === el, model.elements) for el in w.elements]
f(t) = compliance(model, ModelState{eltype(t)}(st.X,
    map(i -> i in idx ? ShellSection(concrete, t[1]) : st.sections[i], eachindex(st.sections))))
Zygote.gradient(f, [0.2])
```

Rigid diaphragms are supported on the pure path. Their master–slave offsets are held at the values of the last `solve!`.
