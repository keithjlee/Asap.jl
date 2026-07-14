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

# global axes
const globalX = SVector(1.0, 0.0, 0.0)
const globalY = SVector(0.0, 1.0, 0.0)
const globalZ = SVector(0.0, 0.0, 1.0)

# ── materials & sections ────────────────────────────────────────────────────
include("Materials/materials.jl")
include("Materials/sections.jl")
export Material, Steel_Nmm, Steel_kNm
export AbstractSection, Section, RigiditySection
export EA, EIx, EIy, GJ, ρA

# ── nodes ───────────────────────────────────────────────────────────────────
include("Nodes/node.jl")
export Node, FIXITIES, fixnode!, planarize!

# ── elements ────────────────────────────────────────────────────────────────
include("Elements/end_conditions.jl")
export EndSprings, EndConditions, rigid_end, pinned_end, RELEASES, release_symbol

include("Elements/kernels/transformation.jl")
include("Elements/kernels/stiffness.jl")
include("Elements/interface.jl")
include("Elements/frame.jl")
include("Elements/variable.jl")
export AbstractElement, FrameElement, TrussElement, VariableElement
export nodes, ndofs, n_internal_dofs, dof_signature, stiffness
export endpoints, midpoint, local_frame
export n_segments, segment_fractions, locate_segment

# ── springs ─────────────────────────────────────────────────────────────────
include("Springs/nodal_springs.jl")
export NodalSpring

# ── loads ───────────────────────────────────────────────────────────────────
include("Loads/loads.jl")
include("Elements/kernels/fixed_end_forces.jl")
export AbstractLoad, NodeLoad, ElementLoad
export NodeForce, NodeMoment, DistributedLoad, LineLoad, TrapezoidLoad, PointLoad, PointMoment, SelfWeight
export fixed_end_forces, condense_fef

# ── model ───────────────────────────────────────────────────────────────────
include("Model/model.jl")
export Model, node_positions, connectivity, volume

# ── analysis ────────────────────────────────────────────────────────────────
include("Analysis/dofs.jl")
include("Analysis/symbolic.jl")
include("Analysis/assemble.jl")
include("Results/results.jl")
include("Analysis/solve.jl")
include("Analysis/functional.jl")
include("Results/recovery.jl")
include("Analysis/cases.jl")
export DofPartition, AnalysisCache
export process!, solve!, assemble_K!, assemble_loads!
export LinearResults, displacement, reaction, element_forces, axial_force
export ElementForceState, LoadTrace, internal_forces, InternalForces
export shear_y, shear_z, torsion, moment_y, moment_z, local_displacements
export LoadCombination, CaseResults, load_cases, solve_cases!, combine, Envelope, envelope
export ModelState, extract_state, solve, compliance, assemble_K

# ── model utilities ─────────────────────────────────────────────────────────
include("Model/utilities.jl")
export clear_supports!, element_connectivity

# ── force density method (self-contained subsystem) ─────────────────────────
include("FDM/FDM.jl")

# ── parametric structure generators (absorbed from AsapToolkit) ──────────────
include("Generation/Generators.jl")
export Frame, SpaceFrame, Warren2D, Pratt2D, SpaceFrameBeam, BakerTruss,
    TrussFrame, GridNetwork, GridFrame
export XGroundStructure, DenseGroundStructure, BoundedGroundStructure
export to_truss, to_frame

# ── plot-ready geometry extraction (absorbed from AsapToolkit) ───────────────
include("Geometry/Displacements.jl")
include("Geometry/Geometry.jl")
export Geo, ModelGeo, TrussGeo, NetworkGeo
export ElementDisplacements, displacements

# ── display (last: defines Base.show for every user-facing type above) ──────
include("ShowMethods.jl")

end # module Asap
