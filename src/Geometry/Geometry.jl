module Geometry

# export functions
export displayGeometry
export meshGeometry
export displayMesh

# export structs
export SubEnclosure
export RayTracingMesh
export GasProperties

# external dependencies
using StaticArrays
using LinearAlgebra
using Plots

# include code
include("SubEnclosure.jl")
include("RayTracingMesh.jl")
include("GasProperties.jl")
include("localWalls.jl")
include("meshGeometry.jl")
include("defineSolidWalls.jl")
include("areaVolumeMesh.jl")
include("displayMesh.jl")
include("displayGeometry.jl")

end