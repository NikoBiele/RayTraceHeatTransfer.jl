module RayTraceHeatTransfer

# external dependencies
using GeometryBasics
using StaticArrays
using LinearAlgebra
using Random
using ProgressMeter
using Measurements
using Base.Threads
using StatsBase
using SparseArrays

# constants
const xVecGlobal2D = SVector(1.0, 0.0)
const yVecGlobal2D = SVector(0.0, 1.0)
const xVecGlobal3D = SVector(1.0, 0.0, 0.0)
const yVecGlobal3D = SVector(0.0, 1.0, 0.0)
const zVecGlobal3D = SVector(0.0, 0.0, 1.0)
const STEFAN_BOLTZMANN = 5.670374419e-8
const h_P = 6.626e-34  # Planck constant
const c0 = 2.998e8     # speed of light  
const k_B = 1.38e-23   # Boltzmann constant
const C2 = h_P * c0 / k_B

# Declare plotting function names so they can be extended by an extension:
function plotMesh end
function plotField end

# throw error if plotting is not enabled
_notenabled(f) = throw(ArgumentError("$f requires the plotting extension. Install and load GLMakie and Plots to enable it,
                                        or ensure that the correct function arguments are provided."))
plotMesh(args...; kwargs...) = _notenabled(:plotMesh)
plotField(args...; kwargs...) = _notenabled(:plotField)

# include code
include(joinpath(@__DIR__, "Domains", "domains.jl"))
include(joinpath(@__DIR__, "Meshing", "meshing.jl"))
include(joinpath(@__DIR__, "RayTracing", "rayTracing.jl"))
include(joinpath(@__DIR__, "HeatTransfer","heatTransfer.jl"))

# export
export PolyVolume2D, PolyFace3D, RayTracingDomain2D, RayTracingDomain3D, ViewFactorDomain3D, viewFactor3D,
        buildSystemMatrices!, solveEquilibrium!, plotMesh, plotField
        
end