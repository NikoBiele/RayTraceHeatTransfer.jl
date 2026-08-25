"""
Regression test: grey vs spectral consistency on a SHARED F matrix.

Pins the receiver-vs-sender b-indexing bug found July 2026 in the grey
solver's absorbed/reflected post-processing (r was computed with b[sender]
instead of b[receiver]; wrong temperatures whenever b is nonuniform across
elements).

The configuration is deliberately the blind spot that hid the bug for years:
NONUNIFORM b (reflective walls eps < 1, non-scattering volumes b = 0) with
participating media. Spectrally uniform properties mean the spectral solve
must reproduce the grey solve EXACTLY -- and by tracing F once and feeding
the identical matrix to both paths, the comparison is deterministic:
machine-precision tolerances, no Monte Carlo allowance.

Any future refactor that reintroduces sender/receiver confusion in either
path fails this test by ~100 K, not by a marginal tolerance.
"""

println("\n" * "-"^60)
println("Testing Grey vs Spectral Consistency (Shared F, Nonuniform b)")
println("-"^60)

using RayTraceHeatTransfer
using Test
using LinearAlgebra
using StatsBase
using StaticArrays
using GeometryBasics
using SparseArrays
using ConvolutionInterpolations

# Deterministic solver-level comparison => tight tolerances.
SHARED_F_TEMP_TOL = 1e-4   # K   (observed agreement ~1e-10; margin for platform variation)
SHARED_F_J_RTOL   = 1e-4   #     (observed ~1e-12)

@testset "Shared-F Grey vs Spectral (nonuniform b)" begin
    T_hot, T_cold  = 1000.0, 0.0
    kappa, sigma_s = 1.0, 0.0
    Ndim, n_bins   = 5, 10
    N_rays         = 1_000_000
    eps_val        = 0.3     # walls b = 0.7, volumes b = 0.0 => nonuniform b
    band_limits    = 10 .^ range(log10(1e-8), log10(0.1), length=n_bins+1)

    vertices   = SVector(Point2(0.0,0.0), Point2(1.0,0.0), Point2(1.0,1.0), Point2(0.0,1.0))
    solidWalls = SVector(true, true, true, true)

    # ---- grey mesh -----------------------------------------------------------
    face_grey = PolyVolume2D{Float64}(vertices, solidWalls, 1, kappa, sigma_s)
    face_grey.T_in_w  = [T_hot, T_cold, T_cold, T_cold]
    face_grey.epsilon = [eps_val, eps_val, eps_val, eps_val]
    face_grey.T_in_g  = -1.0
    mesh_grey = RayTracingDomain2D([face_grey], [(Ndim, Ndim)])

    # ---- spectral mesh, identical physics ------------------------------------
    face_spec = PolyVolume2D{Float64}(vertices, solidWalls, n_bins, kappa, sigma_s)
    face_spec.kappa_g   = fill(kappa, n_bins)
    face_spec.sigma_s_g = fill(sigma_s, n_bins)
    face_spec.q_in_g = 0.0
    face_spec.T_in_g = -1.0
    face_spec.T_in_w = [T_hot, T_cold, T_cold, T_cold]
    face_spec.j_g = zeros(n_bins); face_spec.g_a_g = zeros(n_bins); face_spec.e_g = zeros(n_bins)
    face_spec.r_g = zeros(n_bins); face_spec.g_g   = zeros(n_bins); face_spec.i_g = zeros(n_bins)
    for w in 1:4
        face_spec.j_w[w] = zeros(n_bins); face_spec.g_a_w[w] = zeros(n_bins)
        face_spec.e_w[w] = zeros(n_bins); face_spec.r_w[w]   = zeros(n_bins)
        face_spec.g_w[w] = zeros(n_bins); face_spec.i_w[w]   = zeros(n_bins)
    end
    face_spec.epsilon = [fill(eps_val, n_bins) for _ in 1:4]
    mesh_spec = RayTracingDomain2D([face_spec], [(Ndim, Ndim)])
    mesh_spec.n_spectral_bins        = n_bins
    mesh_spec.wavelength_band_limits = band_limits

    # ---- trace ONCE; both solvers consume the identical F --------------------
    mesh_spec(N_rays; method=:exchange)
    smooth!(mesh_spec)
    F_shared = mesh_spec.F_smooth

    solveEquilibrium!(mesh_grey, F_shared[1])
    solveEquilibrium!(mesh_spec, F_shared)

    # ---- temperature field must agree to solver precision --------------------
    T_grey = [ff.T_g for ff in mesh_grey.fine_mesh[1]]
    T_spec = [ff.T_g for ff in mesh_spec.fine_mesh[1]]

    @test length(T_grey) == length(T_spec)
    @test maximum(abs.(T_spec .- T_grey)) < SHARED_F_TEMP_TOL

    # ---- radiosity (spectral summed over bins) must match grey ---------------
    num_s = length(mesh_spec.surface_mapping)
    num_v = length(mesh_spec.volume_mapping)

    j_grey = zeros(num_s + num_v)
    for ((ci, fi, wi), si) in mesh_grey.surface_mapping
        j_grey[si] = mesh_grey.fine_mesh[ci][fi].j_w[wi]
    end
    for ((ci, fi), vi) in mesh_grey.volume_mapping
        j_grey[num_s + vi] = mesh_grey.fine_mesh[ci][fi].j_g
    end

    j_spec = zeros(num_s + num_v)
    for ((ci, fi, wi), si) in mesh_spec.surface_mapping
        jw = mesh_spec.fine_mesh[ci][fi].j_w[wi]
        @test jw isa AbstractVector           # spectral shape asserted, not assumed
        j_spec[si] = sum(jw)
    end
    for ((ci, fi), vi) in mesh_spec.volume_mapping
        jg = mesh_spec.fine_mesh[ci][fi].j_g
        @test jg isa AbstractVector
        j_spec[num_s + vi] = sum(jg)
    end

    @test all(isapprox.(j_spec, j_grey; rtol=SHARED_F_J_RTOL))

    # ---- the specific quantity the bug corrupted: volume e = j when b = 0 ----
    # (receiver-indexed split => volumes reflect nothing; sender-indexed bug
    #  stripped ~45% here. T^4-consistency pins the post-processing directly.)
    for ((ci, fi), vi) in mesh_spec.volume_mapping
        Tg = mesh_grey.fine_mesh[ci][fi].T_g
        Ts = mesh_spec.fine_mesh[ci][fi].T_g
        if Tg > 1.0
            @test isapprox((Ts / Tg)^4, 1.0; atol=1e-4)
        end
    end
end

println("✓ Shared-F Grey vs Spectral consistency tests complete")