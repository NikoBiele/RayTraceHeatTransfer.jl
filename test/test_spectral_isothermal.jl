# Isothermal enclosure regression tests — spectral solver paths
#
# Physics pinned here: a closed cavity with every wall at a single temperature
# T_wall must solve to T_g == T_wall everywhere, in every spectral mode, for
# any kappa(lambda) and epsilon(lambda). This identity is exact regardless of
# Monte Carlo noise (smoothed F makes uniform T an exact solution), so any
# violation beyond the tolerance is a convention/bookkeeping regression, not
# statistics.
#
# History: v0.11 era bug — the post-convergence temperature writeback in the
# Woodbury solver was handed weighted emission fractions where the e -> T
# inversion requires plain Planck fractions, biasing all radiative-equilibrium
# temperatures by a uniform factor (measured: -2.07 K at 1000 K for
# kappa in [0.5, 1.5] over 5 bins). Fixed by computing plain fractions inside
# the writeback. These tests fail on that bug by three orders of magnitude.
#
# Observed post-fix floors (2e6 rays, 5x5 mesh, 5 bins):
#   uniform path ~2e-13 K, spectral paths ~1e-6 K.
# Tolerance 1e-3 K catches any systematic while staying far above the floors.

using Test
using RayTraceHeatTransfer
using GeometryBasics, StaticArrays

function _isothermal_spectral_error(κ_bins::Vector{Float64},
                                    ε_bins::Vector{Float64};
                                    T_wall = 1000.0, Ndim = 5,
                                    N_rays = 1_000_000,
                                    λ_min = 1e-7, λ_max = 1e-2)
    n_bins = length(κ_bins)
    vertices = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0),
                       Point2(1.0, 1.0), Point2(0.0, 1.0))
    face = PolyVolume2D{Float64}(vertices, SVector(true, true, true, true),
                                 n_bins, 1.0, 0.0)
    face.kappa_g   = copy(κ_bins)
    face.sigma_s_g = fill(0.0, n_bins)
    face.epsilon   = [copy(ε_bins) for _ in 1:4]
    face.T_in_w    = fill(T_wall, 4)
    face.q_in_w    = zeros(4)
    face.T_in_g    = -1.0
    face.q_in_g    = 0.0

    mesh = RayTracingDomain2D([face], [(Ndim, Ndim)])
    mesh.wavelength_band_limits =
        10 .^ range(log10(λ_min), log10(λ_max), length = n_bins + 1)

    mesh(N_rays; method = :exchange)
    smooth!(mesh)
    solveEquilibrium!(mesh, mesh.F_smooth;
                      max_iters = 10_000, convergence_tol = 1e-14)

    return maximum(abs(ff.T_g - T_wall) for ff in mesh.fine_mesh[1])
end

@testset "Isothermal enclosure — spectral solver paths" begin
    tol = 1e-9   # kelvin
    n_bins = 5

    @testset "uniform kappa and epsilon (shared-F direct path)" begin
        err = _isothermal_spectral_error(fill(1.0, n_bins), fill(1.0, n_bins))
        @test err < tol
    end

    @testset "spectral epsilon, uniform kappa (full solver, reflecting)" begin
        # step emissivity: 0.3 below 4 um, 0.9 above (bins 1-3 vs 4-5 at
        # the default log-spaced edges 1e-7..1e-2)
        ε = [0.3, 0.3, 0.3, 0.9, 0.9]
        err = _isothermal_spectral_error(fill(1.0, n_bins), ε)
        @test err < tol
    end

    @testset "spectral kappa, black walls (per-bin F, spectral_variable)" begin
        # this is the configuration that exposed the writeback bug
        κ = collect(range(0.5, 1.5, length = n_bins))
        err = _isothermal_spectral_error(κ, fill(1.0, n_bins))
        @test err < tol
    end
end