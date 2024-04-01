module HeatTransfer

# export user functions
export readMatricesFromCSV
export steadyStateApprox
export steadyStateRigorous
export plotTemperatureField
export validCrosbieSchrenker

# external dependencies
using CSV
using DataFrames
using Plots
using StaticArrays
using LinearAlgebra

# internal dependencies
using ..Geometry
using ..RayTracing

# include code
include("readMatricesFromCSV.jl")
include("steadyStateApprox.jl")
include("steadyStateRigorous.jl")
include("plotTemperatureField.jl")
include("validCrosbieSchrenker.jl")

end