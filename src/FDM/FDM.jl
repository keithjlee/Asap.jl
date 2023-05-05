include("types.jl")
export FDMnode
export FDMelement
export Load
export Network

include("functions.jl")
include("utilities.jl")
export deconstructNode
export process!
export solve!
export forceDensities
export dist
export memberForces
export memberLengths
export FL

export qUpdate!
export xyzUpdate!
export initialLengths
export FDMremoteread
export FDMremotesave