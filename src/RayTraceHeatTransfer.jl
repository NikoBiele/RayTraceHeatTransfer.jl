module RayTraceHeatTransfer

# export functions
export whichEnclosure
export solidWall
export localWalls
export distToSurface
export coefs_FSS_FSG
export coefs_FGS_FGG
export sampleVolume
export lambertSample3D
export sampleSurface
export geometry
export sampleSurfaces
export sampleVolumes
export writeMatricesToCSV
export validCrosbieSchrenker
export steadyStateRigorous
export calculateAreaVolume
export rayTracing_2D
export readMatricesFromCSV

using Plots, LinearAlgebra, CSV, DataFrames, StaticArrays

# include functions
include("whichEnclosure.jl")
include("solidWall.jl")
include("localWalls.jl")
include("distToSurface.jl")
include("coefs_FSS_FSG.jl")
include("coefs_FGS_FGG.jl")
include("sampleVolume.jl")
include("lambertSample3D.jl")
include("sampleSurface.jl")
include("geometry.jl")
include("sampleSurfaces.jl")
include("sampleVolumes.jl")
include("writeMatricesToCSV.jl")
include("validCrosbieSchrenker.jl")
include("steadyStateRigorous.jl")
include("calculateAreaVolume.jl")
include("rayTracing_2D.jl")
include("readMatricesFromCSV.jl")

# example code within this function call
# this function call uses all of the above functions
# except the ones related to CSV-files
# N_rays_tot is the total number of rays traced in the enclosure
N_rays_tot = 10_000_000
# Ndim is the number of divisions of each side of the square enclosure
Ndim = 11
# this function call returns nothing,
# however it prints to the REPL and plots
# it performs a validation of the code against an analytical solution
validCrosbieSchrenker(N_rays_tot,Ndim)
# it is recommenended to save matrices as CSV-files
# especially if they were expensive to obtain
# so that they are not lost if the computer crashes

end
