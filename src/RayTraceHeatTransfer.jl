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
export meshGeometry
export sampleSurfaces
export sampleVolumes
export writeMatricesToCSV
export validCrosbieSchrenker
export steadyStateRigorous
export calculateAreaVolume
export rayTracing_2D
export readMatricesFromCSV
export plotTemperatureField
export sampleDomain
export displayGeometry
# export structs
export TracingMesh
export GasProperties

# include dependencies
using Plots
using LinearAlgebra
using CSV
using DataFrames
using StaticArrays

# include functions
include("Structs.jl")
include("whichEnclosure.jl")
include("solidWall.jl")
include("localWalls.jl")
include("distToSurface.jl")
include("coefs_FSS_FSG.jl")
include("coefs_FGS_FGG.jl")
include("sampleVolume.jl")
include("lambertSample3D.jl")
include("sampleSurface.jl")
include("meshGeometry.jl")
include("sampleSurfaces.jl")
include("sampleVolumes.jl")
include("writeMatricesToCSV.jl")
include("validCrosbieSchrenker.jl")
include("steadyStateRigorous.jl")
include("calculateAreaVolume.jl")
include("rayTracing_2D.jl")
include("readMatricesFromCSV.jl")
include("plotTemperatureField.jl")
include("sampleDomain.jl")
include("displayGeometry.jl")
include("defineSolidWalls.jl")

end
