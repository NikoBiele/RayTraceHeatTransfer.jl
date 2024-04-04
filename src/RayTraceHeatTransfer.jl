module RayTraceHeatTransfer

# export functions and structs for geometry creation
export SubEnclosure # struct
export displayGeometry # function
export meshGeometry # function
export RayTracingMesh # struct
export displayMesh # function
export GasProperties # struct

# export functions for ray tracing
export sampleDomain # function
export writeMatricesToCSV # function

# export functions for heat transfer calculation
export readMatricesFromCSV # function
export steadyState # function (2 methods)
export plotTemperatureField # function
export validCrosbieSchrenker # function

# export viewFactor function
export viewFactor

# include code
include("Geometry/Geometry.jl")
include("RayTracing/RayTracing.jl")
include("HeatTransfer/HeatTransfer.jl")
include("ViewFactor3D/ViewFactor3D.jl")

# internal dependencies
using .Geometry
using .RayTracing
using .HeatTransfer
using .ViewFactor3D

# external dependencies
using StaticArrays
using LinearAlgebra
using Plots
using CSV
using DataFrames

end
