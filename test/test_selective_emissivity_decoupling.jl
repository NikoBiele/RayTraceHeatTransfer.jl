"""
Test: selective-emissivity spectral splitting via exact bin decoupling.

With ALL temperatures prescribed, the K-bin spectral problem decouples
exactly into K independent grey problems: bin k is the grey problem with
wall emissivity eps_k and wall emission eps_k * w_k(T) * A*sigma*T^4, where
w_k is the Planck band fraction. No equilibrium iteration, no coupling.

This test therefore validates the spectral solver's per-bin fields against
K grey solves ON THE SAME shared F matrix, with the grey wall temperature
spoofed so grey emission equals the bin emission:

    T_eff = T_wall * w_k(T_wall)^(1/4),  epsilon = eps_k

Per-bin radiosity must then match to solver precision (not Monte Carlo
tolerance).

This is the first test sensitive to the property-weighted emission
splitting (f_k = eps_k*w_k / sum_i eps_i*w_i), fixed July 2026. The old
plain-Planck splitting (f_k = w_k) produces per-bin errors of order
eps_k/eps_bar - 1 -- tens of percent for the strongly selective epsilon
used here -- so a regression fails loudly, not marginally.

Volumes use bin-uniform kappa so a single traced F serves all bins; the
selectivity lives entirely in epsilon, which is exactly what the weighted
splitting governs.
"""

println("\n" * "-"^60)
println("Testing Selective Emissivity via Exact Bin Decoupling")
println("-"^60)

using RayTraceHeatTransfer
using Test
using LinearAlgebra
using StatsBase
using StaticArrays
using GeometryBasics
using SparseArrays
using ConvolutionInterpolations

DECOUPLE_J_RTOL = 1e-4    # deterministic shared-F comparison; observed ~1e-10

@testset "Selective epsilon: per-bin vs independent grey solves" begin
    # ---- configuration -------------------------------------------------------
    T_walls = [1000.0, 300.0, 300.0, 300.0]
    T_gas   = 600.0
    kappa, sigma_s = 1.0, 0.0
    Ndim    = 4
    n_bins  = 3
    N_rays  = 1_000_000
    # Strongly selective: old plain-w_k splitting errs by ~eps_k/eps_bar - 1
    eps_bins = [0.9, 0.5, 0.2]
    # Band edges bracketing the Planck peaks of 300-1000 K (peak lambda*T ~ 2.9e-3 m*K)
    band_limits = [1.0e-6, 3.0e-6, 8.0e-6, 1.0e-3]

    vertices   = SVector(Point2(0.0,0.0), Point2(1.0,0.0), Point2(1.0,1.0), Point2(0.0,1.0))
    solidWalls = SVector(true, true, true, true)

    # ---- spectral mesh: selective walls, uniform kappa, ALL T prescribed -----
    face_spec = PolyVolume2D{Float64}(vertices, solidWalls, n_bins, kappa, sigma_s)
    face_spec.kappa_g   = fill(kappa, n_bins)
    face_spec.sigma_s_g = fill(sigma_s, n_bins)
    face_spec.q_in_g = 0.0
    face_spec.T_in_g = T_gas                     # prescribed => exact decoupling
    face_spec.T_in_w = copy(T_walls)
    face_spec.j_g = zeros(n_bins); face_spec.g_a_g = zeros(n_bins); face_spec.e_g = zeros(n_bins)
    face_spec.r_g = zeros(n_bins); face_spec.g_g   = zeros(n_bins); face_spec.i_g = zeros(n_bins)
    for w in 1:4
        face_spec.j_w[w] = zeros(n_bins); face_spec.g_a_w[w] = zeros(n_bins)
        face_spec.e_w[w] = zeros(n_bins); face_spec.r_w[w]   = zeros(n_bins)
        face_spec.g_w[w] = zeros(n_bins); face_spec.i_w[w]   = zeros(n_bins)
    end
    face_spec.epsilon = [copy(eps_bins) for _ in 1:4]
    mesh_spec = RayTracingDomain2D([face_spec], [(Ndim, Ndim)])
    mesh_spec.n_spectral_bins        = n_bins
    mesh_spec.wavelength_band_limits = band_limits
    # spectral_mode left to the package's own detection where possible; set if required:
    @test mesh_spec.spectral_mode == :spectral_variable

    # ---- trace ONCE ----------------------------------------------------------
    mesh_spec(N_rays; method=:exchange)
    F_shared = mesh_spec.F_smooth

    solveEquilibrium!(mesh_spec, F_shared)

    num_s = length(mesh_spec.surface_mapping)
    num_v = length(mesh_spec.volume_mapping)
    n_elem = num_s + num_v

    # collect spectral per-bin radiosity, element-ordered
    j_spec_bins = zeros(n_elem, n_bins)
    for ((ci, fi, wi), si) in mesh_spec.surface_mapping
        jw = mesh_spec.fine_mesh[ci][fi].j_w[wi]
        @test jw isa AbstractVector
        @test length(jw) == n_bins
        j_spec_bins[si, :] .= jw
    end
    for ((ci, fi), vi) in mesh_spec.volume_mapping
        jg = mesh_spec.fine_mesh[ci][fi].j_g
        @test jg isa AbstractVector
        j_spec_bins[num_s + vi, :] .= jg
    end

    # ---- K independent grey reference solves on the SAME F -------------------
    # Band fractions from the package itself (same C2, same normalization as
    # the solver) so the reference is internally consistent at tight tolerance.
    w_at = Dict{Float64, Vector{Float64}}()
    for T in unique(vcat(T_walls, T_gas))
        w_at[T] = RayTraceHeatTransfer.getBinsEmissionFractions(mesh_spec, fill(T, n_elem))[1, :]
    end

    for bin in 1:n_bins
        T_w_eff = [Tw * w_at[Tw][bin]^0.25 for Tw in T_walls]
        T_g_eff = T_gas * w_at[T_gas][bin]^0.25

        face_g = PolyVolume2D{Float64}(vertices, solidWalls, 1, kappa, sigma_s)
        face_g.T_in_w  = T_w_eff
        face_g.epsilon = fill(eps_bins[bin], 4)
        face_g.T_in_g  = T_g_eff
        mesh_g = RayTracingDomain2D([face_g], [(Ndim, Ndim)])

        solveEquilibrium!(mesh_g, F_shared[1])     # no trace: shared F

        j_grey = zeros(n_elem)
        for ((ci, fi, wi), si) in mesh_g.surface_mapping
            j_grey[si] = mesh_g.fine_mesh[ci][fi].j_w[wi]
        end
        for ((ci, fi), vi) in mesh_g.volume_mapping
            j_grey[num_s + vi] = mesh_g.fine_mesh[ci][fi].j_g
        end

        @testset "bin $bin (eps = $(eps_bins[bin]))" begin
            @test all(isapprox.(j_spec_bins[:, bin], j_grey;
                                rtol=DECOUPLE_J_RTOL, atol=1e-4))
        end
    end
end

println("✓ Selective Emissivity Decoupling tests complete")