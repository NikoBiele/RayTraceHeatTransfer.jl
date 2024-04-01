module ViewFactor3D

# export functions
export viewFactor

# external dependencies
using LinearAlgebra

# include code
include("viewFactor.jl")
include("Cl.jl")
include("imagLi_2.jl")
include("edgePairParameters.jl")
include("fparallel.jl")
include("f.jl")

end