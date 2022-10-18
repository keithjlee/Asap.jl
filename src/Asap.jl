module Asap

using Base: StackFrame, Float64, Symbol
using Statistics
using LinearAlgebra
using SparseArrays
using GeometryBasics
using FileIO
using ColorSchemes
using Colors
using GLMakie

# data types used for analysis
include("types.jl")
export Node
export Element
export Load
export Structure
export Geometry

# structural analysis
include("R.jl") # rotation matrices
include("K.jl") # stiffness matrices
include("structuralAnalysis.jl") #operations
export analyze!
export dofMaker
export pseudoSize
export pseudoSize!
export structureMass

# method extensions/general QOL utilities
using JSON
include("utilities.jl")
export karamba2asap

# viz tools
include("visualization.jl")
export axo
export structurePlot

end