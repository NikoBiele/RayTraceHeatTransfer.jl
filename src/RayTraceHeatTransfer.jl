module RayTraceHeatTransfer

const STEFAN_BOLTZMANN = 5.670374419e-8

# external dependencies
using GLMakie
using GeometryBasics
using StaticArrays
using LinearAlgebra
using Plots
using Random
using ProgressMeter
using Measurements
using Base.Threads
using StatsBase

# export
# 2D participating media
export PolyFace2D, RayTracingMeshOptim, plotMesh2D,
        GasProperties, steadyState2D!, plotField2D,
# 3D transparent medium
        Domain3D_faces, plotMesh3D, steadyState3D!,
        plotField3D

# include code
include("./Domains/domains.jl")
include("./HeatTransfer/heatTransfer.jl")
include("./RayTracing/rayTracing.jl")

end