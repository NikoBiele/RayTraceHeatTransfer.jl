module RayTraceHeatTransfer

const STEFAN_BOLTZMANN = 5.670374419e-8
const h_P = 6.626e-34  # Planck constant
const c0 = 2.998e8     # speed of light  
const k_B = 1.38e-23   # Boltzmann constant
const C2 = h_P * c0 / k_B

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

# Declare plotting function names so they can be extended by an extension:
function plotMesh2D end
function plotField2D end
function plotMesh3D end
function plotField3D end

# throw error if plotting is not enabled
_notenabled(f) = throw(ArgumentError("$f requires the plotting extension. Install and load GLMakie and Plots to enable it,
                                        or ensure that the correct function arguments are provided."))
plotMesh2D(args...; kwargs...) = _notenabled(:plotMesh2D)
plotField2D(args...; kwargs...) = _notenabled(:plotField2D)
plotMesh3D(args...; kwargs...) = _notenabled(:plotMesh3D)
plotField3D(args...; kwargs...) = _notenabled(:plotField3D)

# include code
include(joinpath(@__DIR__, "Domains",    "domains.jl"))
include(joinpath(@__DIR__, "HeatTransfer","heatTransfer.jl"))
include(joinpath(@__DIR__, "RayTracing", "rayTracing.jl"))

# export
# 2D participating media
export PolyFace2D, RayTracingDomain2D, plotMesh2D, build_system_matrices2D!,
        steadyStateGrey2D!, steadyStateSpectral2D!, plotField2D,
# 3D transparent medium
        ViewFactorDomain3D, plotMesh3D, viewFactor, build_system_matrices3D!,
        steadyStateGrey3D!, steadyStateSpectral3D!,
        plotField3D
        
end