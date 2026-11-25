"""
    Asap

Structural analysis for trusses and frames — the v1.0 core.

Three separated layers: definition (`Model`, `Node`, elements, loads,
springs), analysis structure (`AnalysisCache`, built by `process!`), and
results (`LinearResults`, returned by `solve!`). One set of pure element
kernels feeds both an in-place zero-allocation assembly path and a pure
functional path (`ModelState` → `solve`) that automatic differentiation
engines traverse natively (see `ext/AsapChainRulesExt.jl`).

See `docs/MODERNIZATION.md` for the architecture and roadmap.
"""
module Asap

using LinearAlgebra, SparseArrays, StaticArrays

# solve/solve! extend the CommonSolve verbs, so loading other CommonSolve
# consumers (LinearSolve, the SciML stack) NEVER shadows Asap's API — both
# packages export the same binding
import CommonSolve: solve, solve!

# global axes
const globalX = SVector(1.0, 0.0, 0.0)
const globalY = SVector(0.0, 1.0, 0.0)
const globalZ = SVector(0.0, 0.0, 1.0)

# ── materials & sections ────────────────────────────────────────────────────
include("Materials/materials.jl")
include("Materials/sections.jl")
include("Materials/shell_sections.jl")
export Material, Steel_Nmm, Steel_kNm
export AbstractSection, Section, RigiditySection
export EA, EIx, EIy, GJ, ρA
export AbstractShellSection, ShellSection, ShellRigiditySection
export membrane_rigidity, bending_rigidity, shear_rigidity, drilling_rigidity, areal_mass, thickness

# ── nodes ───────────────────────────────────────────────────────────────────
include("Nodes/node.jl")
export Node, FIXITIES, fixnode!, planarize!

# ── elements ────────────────────────────────────────────────────────────────
include("Elements/end_conditions.jl")
export EndSprings, EndConditions, rigid_end, pinned_end, RELEASES, release_symbol

include("Elements/kernels/transformation.jl")
include("Elements/kernels/stiffness.jl")
include("Elements/kernels/shell.jl")
include("Elements/kernels/shell_tri.jl")
include("Elements/interface.jl")
include("Elements/frame.jl")
include("Elements/variable.jl")
include("Elements/shell.jl")
include("Elements/link.jl")
export AbstractElement, AbstractLineElement, FrameElement, TrussElement, VariableElement
export nodes, ndofs, n_internal_dofs, dof_signature, stiffness
export endpoints, midpoint, local_frame
export n_segments, segment_fractions, locate_segment
export ShellElement, shell_frame, shell_stiffness, itw_membrane_stiffness
export andes_membrane_stiffness, andes_beta0, mitc3p_plate_stiffness
export shell_resultants, shell_envelope, section_cut, LinkElement

# ── springs ─────────────────────────────────────────────────────────────────
include("Springs/nodal_springs.jl")
export NodalSpring

# ── loads ───────────────────────────────────────────────────────────────────
include("Loads/loads.jl")
include("Elements/kernels/fixed_end_forces.jl")
include("Loads/shell_loads.jl")
include("Elements/kernels/shell_loads_recovery.jl")
export AbstractLoad, NodeLoad, ElementLoad
export NodeForce, NodeMoment, DistributedLoad, LineLoad, TrapezoidLoad, PointLoad, PointMoment, SelfWeight
export fixed_end_forces, condense_fef
export AreaLoad

# ── model ───────────────────────────────────────────────────────────────────
include("Model/diaphragm.jl")
include("Model/model.jl")
export Model, node_positions, connectivity, volume
export RigidDiaphragm

# ── analysis ────────────────────────────────────────────────────────────────
include("Analysis/dofs.jl")
include("Analysis/constraints.jl")   # ConstraintMap is a field type of AnalysisCache
include("Analysis/symbolic.jl")
include("Analysis/assemble.jl")
include("Results/results.jl")
include("Analysis/solve.jl")
include("Analysis/functional.jl")
include("Results/recovery.jl")
include("Analysis/cases.jl")
include("Results/shell_results.jl")    # after cases.jl: shell_envelope takes CaseResults
export DofPartition, AnalysisCache
export process!, solve!, assemble_K!, assemble_loads!
export CachedSolver
export LinearResults, displacement, reaction, element_forces, axial_force
export ElementForceState, LoadTrace, internal_forces, InternalForces
export shear_y, shear_z, torsion, moment_y, moment_z, local_displacements
export LoadCombination, CaseResults, load_cases, solve_cases!, combine, Envelope, envelope
export ModelState, extract_state, solve, compliance, assemble_K, assemble_Pf, fef_local
export condensed_fef, element_fef_global

# ── model utilities ─────────────────────────────────────────────────────────
include("Model/utilities.jl")
export clear_supports!, element_connectivity

# ── force density method (self-contained subsystem) ─────────────────────────
include("FDM/FDM.jl")

# ── parametric structure generators (absorbed from AsapToolkit) ──────────────
include("Generation/Generators.jl")
include("Generation/ShellMesh.jl")
include("Generation/Gmsh.jl")
include("Generation/WallMesh.jl")
export Frame, SpaceFrame, Warren2D, Pratt2D, SpaceFrameBeam, BakerTruss,
    TrussFrame, GridNetwork, GridFrame
export WallMesh, mesh_wall, wall_cut, edge_loads
export AbstractShellMesh, ShellMesh, shell_mesh, mesh_polygon, boundary_nodes, nodes_on_segment,
    elements_in, total_reaction, circle_polygon, group_nodes, group_elements, read_gmsh
export XGroundStructure, DenseGroundStructure, BoundedGroundStructure
export to_truss, to_frame

# ── plot-ready geometry extraction (absorbed from AsapToolkit) ───────────────
include("Geometry/Displacements.jl")
include("Geometry/Geometry.jl")
include("Geometry/ShellGeometry.jl")
export Geo, ModelGeo, TrussGeo, NetworkGeo, ShellGeo, nodal_average, principal_forces, shell_faces, shell_face_owner
export ElementDisplacements, displacements

# ── display (last: defines Base.show for every user-facing type above) ──────
include("ShowMethods.jl")

end # module Asap
