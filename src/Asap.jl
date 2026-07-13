module Asap

using LinearAlgebra, SparseArrays

# global axes
const globalX::Vector{Float64} = [1., 0., 0.]
const globalY::Vector{Float64} = [0., 1., 0.]
const globalZ::Vector{Float64} = [0., 0., 1.]

include("legacy/Materials_Sections/material.jl")
include("legacy/Materials_Sections/section.jl")
export Material
export Section
export TrussSection
export Steel_Nmm
export Steel_kNm

include("legacy/Nodes/nodes.jl")
include("legacy/Nodes/utilities.jl")
export Node
export TrussNode
export planarize!
export fixnode!

include("legacy/Elements/elements.jl")
include("legacy/Elements/K.jl")
include("legacy/Elements/R.jl")
include("legacy/Elements/utilities.jl")
export Element
export BridgeElement
export TrussElement
export release!
export endpoints
export midpoint
export axial_force

include("legacy/Loads/loads.jl")
include("legacy/Loads/utilities.jl")
include("legacy/Loads/fixedEndForces.jl")
export NodeForce
export NodeMoment
export LineLoad
export GravityLoad
export PointLoad

include("legacy/Model/model.jl")
include("legacy/Model/utilities.jl")
export Model
export TrussModel
export update_DOF!
export connectivity
export node_positions
export volume


include("legacy/Model/preprocessing.jl")
include("legacy/Model/bridgeprocessing.jl")
include("legacy/Model/postprocessing.jl")
include("legacy/Model/analysis.jl")
export process!
export solve!
export solve
export connectivity
export node_positions
export volume

# FORCE DENSITY METHOD
include("FDM/FDM.jl")


end 