module Asap

using LinearAlgebra, SparseArrays, Statistics

# global axes
const globalX::Vector{Float64} = [1., 0., 0.]
const globalY::Vector{Float64} = [0., 1., 0.]
const globalZ::Vector{Float64} = [0., 0., 1.]

include("Materials_Sections/material.jl")
include("Materials_Sections/section.jl")
export Material
export Section
export TrussSection
export Steel_Nmm
export Steel_kNm

include("Nodes/nodes.jl")
include("Nodes/utilities.jl")
export Node
export TrussNode
export planarize!
export fixnode!

include("Elements/elements.jl")
include("Elements/K.jl")
include("Elements/R.jl")
include("Elements/utilities.jl")
export Element
export BridgeElement
export TrussElement
export release!
export endpoints
export midpoint
export axial_force

include("Loads/loads.jl")
include("Loads/utilities.jl")
include("Loads/fixedEndForces.jl")
export NodeForce
export NodeMoment
export LineLoad
export GravityLoad
export PointLoad

include("Model/model.jl")
include("Model/utilities.jl")
export Model
export TrussModel
export update_DOF!
export connectivity
export node_positions
export volume


include("Model/preprocessing.jl")
include("Model/bridgeprocessing.jl")
include("Model/postprocessing.jl")
include("Model/analysis.jl")
export process!
export solve!
export solve
export connectivity

# FORCE DENSITY METHOD
include("FDM/FDM.jl")


end 