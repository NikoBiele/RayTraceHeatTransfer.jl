module RayTraceHeatTransfer

# external dependencies
using GeometryBasics
using StaticArrays
using LinearAlgebra
using Random
using ProgressMeter
using Base.Threads
using StatsBase
using SparseArrays
using Krylov
using ConvolutionInterpolations
using Sobol: Sobol, SobolSeq

# constants
const xVecGlobal2D = SVector(1.0, 0.0)
const yVecGlobal2D = SVector(0.0, 1.0)
const xVecGlobal3D = SVector(1.0, 0.0, 0.0)
const yVecGlobal3D = SVector(0.0, 1.0, 0.0)
const zVecGlobal3D = SVector(0.0, 0.0, 1.0)

# Exact SI defining constants (2019 redefinition of the SI)
const h_P = 6.62607015e-34   # Planck constant [J s]
const c0  = 299792458.0      # speed of light in vacuum [m/s]
const k_B = 1.380649e-23     # Boltzmann constant [J/K]

"""
    stefan_boltzmann(T)

Stefan–Boltzmann constant σ = 2π⁵k⁴/(15h³c²) in the floating-point type `T`.
Since 2019, h, c and k are exact by definition, so σ is known to any precision.
"""
stefan_boltzmann(::Type{Float64}) = 5.6703744191844294e-8      # correctly rounded
stefan_boltzmann(::Type{T}) where {T} =
    2 * T(pi)^5 * (T(1380649) / T(10)^29)^4 /                   # 2π⁵k⁴
    (15 * (T(662607015) / T(10)^42)^3 * T(299792458)^2)         # 15h³c²

const STEFAN_BOLTZMANN = stefan_boltzmann(Float64)  # [W/(m² K⁴)]
const C1 = 1.1910429723971884e-16   # 2 h c² [W m² / sr]
const C2 = 1.4387768775039337e-2    # h c / k [m K]

# Declare plotting function names so they can be extended by an extension:
function plotMesh end
function plotField end
function plotDomain3D end

# throw error if plotting is not enabled
_notenabled(f) = throw(ArgumentError("$f requires the plotting extension. Install and load GLMakie and Plots to enable it,
                                        or ensure that the correct function arguments are provided."))
plotMesh(args...; kwargs...) = _notenabled(:plotMesh)
plotField(args...; kwargs...) = _notenabled(:plotField)
plotDomain3D(args...; kwargs...) = _notenabled(:plotDomain3D)

# include code
include(joinpath(@__DIR__, "Domains", "domains.jl"))
include(joinpath(@__DIR__, "Meshing", "meshing.jl"))
include(joinpath(@__DIR__, "RayTracing", "rayTracing.jl"))
include(joinpath(@__DIR__, "Smoothing", "smooth.jl"))
include(joinpath(@__DIR__, "HeatTransfer","heatTransfer.jl"))

# export
export PolyVolume2D,
       PolyFace3D, 
       RayTracingDomain2D,
       ViewFactorDomain3D,
       viewFactor3D,
       buildSystemMatrix,
       solveEquilibrium!,
       plotMesh,
       plotDomain3D,
       plotField,
       RayRecorder,
       collect_rays,
       smooth!,
       RayTracingDomain3D_surfaces,
       elementFaces,
       elementParents,
       AbstractSpectralModel,
       PlanckBands,
       ConstantWeights,
       PiecewiseBands,
       adaptiveSpectralBins,
       exchangeFactors!,
       AbstractDirectionalModel,
       AngularBins,
       IsotropicScattering,
       HenyeyGreenstein,
       TabulatedScattering,
       DiffuseReflection,
       SpecularReflection,
       TabulatedReflection
end