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
export areaVolumeMesh
export rayTracing_CPU
export readMatricesFromCSV
export plotTemperatureField
export sampleDomain
export displayMesh
# export structs
export TracingMesh
export GasProperties

# include dependencies
using Plots
using LinearAlgebra
using CSV
using DataFrames
using StaticArrays

# include code
include("structs.jl")
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
include("areaVolumeMesh.jl")
include("rayTracing_CPU.jl")
include("readMatricesFromCSV.jl")
include("plotTemperatureField.jl")
include("sampleDomain.jl")
include("displayMesh.jl")
include("defineSolidWalls.jl")

end
