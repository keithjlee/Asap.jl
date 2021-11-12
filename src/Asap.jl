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

include("types.jl")
export Node
export Element
export Load
export Structure
export Geometry

include("structuralAnalysis.jl")
export analyze
export dofMaker
export pseudoSize
export pseudoSize!
export structureMass

include("visualization.jl")
export axo
export structurePlot

end