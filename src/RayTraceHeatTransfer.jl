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
export plotTrapezoids
export sampleDomain

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
include("plotTrapezoids.jl")
include("sampleDomain.jl")

end
