module RayTracing

# import dependencies
using Plots

# export functions
export sampleDomain
export writeMatricesToCSV

# external dependencies
using StaticArrays
using LinearAlgebra
using CSV
using DataFrames

# internal dependencies
using ..Geometry

# include code
include("coefs_exchange.jl")
include("lambertSample3D.jl")
include("rayTracing_CPU.jl")
include("sampleDomain.jl")
include("sampleSurface.jl")
include("sampleSurfaces.jl")
include("sampleVolume.jl")
include("sampleVolumes.jl")
include("solidWalls.jl")
include("distToSurface.jl")
include("whichCell.jl")
include("whichSubEnclosure.jl")
include("isotropicScatter.jl")
include("writeMatricesToCSV.jl")

end