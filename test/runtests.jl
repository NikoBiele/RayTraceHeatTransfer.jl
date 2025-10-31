using RayTraceHeatTransfer
using Test
using LinearAlgebra
using StatsBase
using StaticArrays
using GeometryBasics

const VF_TOLERANCE = 1e-5 # Tolerance for view factor comparisons
const TEMP_TOLERANCE = 5.0  # K, Tolerance for absolute temperature
const ENERGY_TOLERANCE = 1e-4 # W, Absolute tolerance for energy balance
const ANALYTICAL_TOLERANCE = 0.05 # 5% tolerance vs analytical solution (to avoid too much sampling)
const SPECTRAL_TOLERANCE = 0.05 # 5% tolerance for spectral comparisons
const CONSISTENCY_TOLERANCE = 0.05  # 5% tolerance for consistency checks (ray tracing vs. exchange factors)

@testset "RayTraceHeatTransfer.jl" begin
    println("\n" * "="^80)
    println("STARTING COMPREHENSIVE TEST SUITE")
    println("="^80)
    
    # Test 3D transparent surfaces (view factors and heat transfer)
    @testset "3D Surfaces (Transparent)" begin
        include("test_3d_viewfactors.jl")
        include("test_3d_heat_transfer.jl")
    end
    
    # Test 2D participating media
    @testset "2D Participating Media" begin
        include("test_2d_grey.jl")
        include("test_2d_spectral.jl")
    end
    
    # Test spectral capabilities
    @testset "Spectral Radiation" begin
        include("test_spectral_consistency.jl")
    end
    
    println("\n" * "="^80)
    println("TEST SUITE COMPLETE")
    println("="^80)
end