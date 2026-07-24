"""
Test 2D grey radiative transfer with NON-ZERO surface reflection.

The other grey tests use fully absorbing surfaces (epsilon = 1), which
makes B = 0 and bypasses the full A/R/C/D matrix structure of the
GERT method. This file exercises the reflection-scattering path:

- Energy conservation with reflecting walls and a gas in radiative
  equilibrium (the configuration that distinguishes the corrected
  method from a buggy one)
- Parallel-plate-approximation against the textbook closed-form for
  two infinite gray-diffuse plates
"""

println("\n" * "-"^60)
println("Testing 2D Grey with Surface Reflection")
println("-"^60)

using RayTraceHeatTransfer
using Test
using LinearAlgebra
using StatsBase
using StaticArrays
using GeometryBasics
using SparseArrays
using ConvolutionInterpolations

ENERGY_TOLERANCE = 1e-4 # W, Absolute tolerance for energy balance
ANALYTICAL_TOLERANCE = 0.05 # 5% tolerance vs analytical solution (to avoid too much sampling)
const SIGMA = 5.670374419e-8  # Stefan-Boltzmann constant, W/(m^2 K^4)

#############################################################################
### TEST 1: ENERGY CONSERVATION WITH REFLECTING WALLS ######################
#############################################################################
#
# Unit square with one hot wall (bottom, T = 1000 K) and three cold walls
# (T = -1, q = 0) with reflectivity rho = 0.5 (i.e., epsilon = 0.5). Gas
# in radiative equilibrium with kappa = 1, sigma_s = 0. Total energy in
# must equal total energy out at steady state, regardless of the
# reflectivity distribution.

@testset "Energy Conservation with Reflecting Walls" begin
    T_hot  = 1000.0
    Ndim   = 11
    N_rays_total = 1_000_000

    vertices = SVector(
        Point2(0.0, 0.0),
        Point2(1.0, 0.0),
        Point2(1.0, 1.0),
        Point2(0.0, 1.0)
    )
    solidWalls = SVector(true, true, true, true)
    face = PolyVolume2D{Float64}(vertices, solidWalls, 1, 1.0, 0.0)

    # Bottom wall hot (prescribed T), other walls reflecting and in
    # radiative equilibrium (prescribed q = 0).
    face.T_in_w  = [T_hot, -1.0, -1.0, -1.0]
    face.q_in_w  = [0.0,    0.0,  0.0,  0.0]
    face.epsilon = [1.0,    0.5,  0.5,  0.5]
    face.T_in_g  = -1.0
    face.q_in_g  = 0.0

    mesh = RayTracingDomain2D([face], [(Ndim, Ndim)])
    mesh(N_rays_total; method = :exchange)
    solveEquilibrium!(mesh, mesh.F_smooth)

    @test abs(mesh.energy_error) < ENERGY_TOLERANCE
end

#############################################################################
### TEST 2: PARALLEL-PLATE APPROXIMATION AGAINST TEXTBOOK ##################
#############################################################################
#
# Wide thin enclosure (width 100, height 1) approximates two infinite
# parallel plates in the central region (far from the narrow sidewalls).
# Hot wall at the bottom, cold wall at the top, both gray-diffuse with
# emissivity 0.5. Sidewalls fully absorbing at T = 0. Transparent medium.
#
# Textbook: net flux from hot plate to cold plate per unit area is
#     q_12 = sigma * (T_hot^4 - T_cold^4) / (1/eps + 1/eps - 1)
# At central hot-wall elements (far from sidewalls), the computed source
# flux per area should match this within ray-tracing accuracy.

@testset "Parallel Plates vs Textbook" begin
    eps_plates = 0.5
    T_hot  = 1000.0
    T_cold =    0.0
    W = 100.0
    H = 1.0
    Nx = 21
    Ny = 2
    N_rays_total = 10_000_000

    vertices = SVector(
        Point2(0.0, 0.0),
        Point2(W,   0.0),
        Point2(W,   H),
        Point2(0.0, H)
    )
    solidWalls = SVector(true, true, true, true)
    face = PolyVolume2D{Float64}(vertices, solidWalls, 1, 1e-3, 0.0) # almost transparent (to help smoothing)

    # Wall numbering: 1 = bottom (hot), 2 = right, 3 = top (cold), 4 = left
    face.T_in_w  = [T_hot,      T_cold, T_cold,     T_cold]
    face.q_in_w  = [0.0,        0.0,    0.0,        0.0]
    face.epsilon = [eps_plates, 1.0,    eps_plates, 1.0]
    face.T_in_g  = -1.0
    face.q_in_g  = 0.0

    mesh = RayTracingDomain2D([face], [(Nx, Ny)])
    mesh(N_rays_total; method = :exchange, k_dykstra=500)
    solveEquilibrium!(mesh, mesh.F_smooth)

    @test abs(mesh.energy_error) < ENERGY_TOLERANCE

    # Textbook prediction (per unit area)
    q_textbook = SIGMA * (T_hot^4 - T_cold^4) / (1/eps_plates + 1/eps_plates - 1)

    # Extract central hot-wall elements (the bottom row of fine faces, in
    # the middle of the width, where the parallel-plate approximation
    # holds). The fine-mesh ordering for a single coarse face with
    # (Nx, Ny) subdivision is column-major: bottom-row hot-wall elements
    # are fine indices 1 through Nx, wall_idx = 1.
    n_central = max(1, Nx ÷ 5)
    center_lo = (Nx - n_central) ÷ 2 + 1
    center_hi = center_lo + n_central - 1

    q_per_area_samples = Float64[]
    for col in center_lo:center_hi
        f = mesh.fine_mesh[1][col]
        push!(q_per_area_samples, f.q_w[1] / f.area[1])
    end
    q_per_area_mean = sum(q_per_area_samples) / length(q_per_area_samples)

    @test abs(q_per_area_mean - q_textbook) / q_textbook < ANALYTICAL_TOLERANCE
end

println("✓ 2D Grey with Surface Reflection tests complete")