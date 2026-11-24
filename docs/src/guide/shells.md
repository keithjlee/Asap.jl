# Shells, walls and diaphragms

```@setup shells
using Asap
```

Asap's [`ShellElement`](@ref) is a 4-node flat shell: a drilling-rotation membrane (Ibrahimbegović–Taylor–Wilson, with two complementary incompatible modes) combined with the MITC4 Reissner–Mindlin plate. It mixes freely with frames and trusses, and its drilling DOF lets an in-plane beam transfer moment into a wall at a shared node.

The formulation and its verification (patch tests, MacNeal–Harder, Timoshenko & Woinowsky-Krieger plates, the Scordelis–Lo roof) are documented in `docs/shells/DERIVATION.md` in the repository.

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
- [`nodal_average`](@ref) smooths any per-element field to the nodes for contour plots.
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

For fast reverse mode, load ForwardDiff alongside Zygote. The shell kernels' rules (`AsapShellADExt`) evaluate each element's local Jacobian with ForwardDiff instead of tracing it, about 15× faster. Mooncake (`AsapShellMooncakeExt`) and Enzyme reverse mode (`AsapShellEnzymeExt`) use the same rules; both fail on the shell kernels natively. All backends agree to 1e-16; Enzyme is the fastest (see `examples/ad/ad_backends_shells.jl`).

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
