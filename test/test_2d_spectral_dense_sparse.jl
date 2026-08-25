println("\n" * "-"^60)
println("Testing 2D spectral variable (dense/sparse Woodbury)")
println("-"^60)

using RayTraceHeatTransfer
using Test
using LinearAlgebra
using StatsBase
using StaticArrays
using GeometryBasics
using SparseArrays
using ConvolutionInterpolations

const T_HOT_CS = 1000.0

# ---- Crosbie & Schrenker (1984) tabulated reference --------------------------
const TAU_REF_CS = [0.0, 0.00611, 0.02037, 0.04251, 0.07216, 0.10884, 0.15194,
                    0.20076, 0.25449, 0.31225, 0.37309, 0.43602, 0.50000, 0.56398,
                    0.62691, 0.68775, 0.74551, 0.79924, 0.84806, 0.89116, 0.92784,
                    0.95749, 0.97963, 0.99390, 1.00000]

const S_REF_CS = [0.6293, 0.6198, 0.6017, 0.5767, 0.5460, 0.5108, 0.4724,
                  0.4323, 0.3919, 0.3525, 0.3153, 0.2810, 0.2500, 0.2224,
                  0.1981, 0.1768, 0.1584, 0.1424, 0.1287, 0.1171, 0.1073,
                  0.0992, 0.0930, 0.0885, 0.0863]

const S_INTERP_CS = convolution_interpolation((TAU_REF_CS,), S_REF_CS; kernel = :b7)

centerline_taus_cs(N_side) = collect(range(1 / (2N_side), 1 - 1 / (2N_side), length = N_side))

# ---- problem builder ---------------------------------------------------------
function build_cs_spectral(N_side::Int, n_bins::Int)
    vertices = SVector(
        Point2(0.0, 0.0), Point2(1.0, 0.0),
        Point2(1.0, 1.0), Point2(0.0, 1.0))
    solidWalls = SVector(true, true, true, true)
    face = PolyVolume2D{Float64}(vertices, solidWalls, n_bins, 1.0, 0.0)
    # ~1% monotonic kappa perturbation: triggers :spectral_variable (Woodbury)
    face.kappa_g = [1.0 * (1.0 + 0.01 * (k - 1) / (n_bins - 1)) for k in 1:n_bins]
    face.sigma_s_g = zeros(n_bins)
    face.T_in_w  = [T_HOT_CS, 0.0, 0.0, 0.0]
    face.epsilon = [fill(1.0, n_bins) for _ in 1:4]
    face.T_in_g  = -1.0
    face.q_in_g  = 0.0
    mesh = RayTracingDomain2D([face], [(N_side, N_side)])
    mesh.wavelength_band_limits = 10 .^ range(log10(1e-8), log10(0.1), length = n_bins + 1)
    return mesh
end

# ---- extraction + error ------------------------------------------------------
function centerline_rms_cs(mesh, N_side)
    all_T = [ff.T_g for ff in mesh.fine_mesh[1]]
    T_center = reshape(all_T, N_side, N_side)[(N_side - 1) ÷ 2 + 1, :]
    S_comp = (T_center ./ T_HOT_CS) .^ 4
    S_ref  = [S_INTERP_CS(t) for t in centerline_taus_cs(N_side)]
    return sqrt(mean((S_comp .- S_ref) .^ 2))
end

# =============================================================================
@testset "Crosbie-Schrenker spectral tau=1 (dense/sparse Woodbury)" begin

    N_side, n_bins = 11, 10
    N_rays = (4 * N_side + N_side^2) * 1_000

    @testset "Sparse/dense spectral C&S" begin
        mesh = build_cs_spectral(N_side, n_bins)
        mesh(N_rays; method = :exchange)
        smooth!(mesh)

        @test mesh.spectral_mode == :spectral_variable        # Woodbury dispatch
        @test all(count(<(0.0), F) == 0 for F in mesh.F_smooth)

        # raw solve (sparse)
        solveEquilibrium!(mesh, mesh.F_raw)
        err_raw = centerline_rms_cs(mesh, N_side)

        # smoothed solve (sparse or dense)
        solveEquilibrium!(mesh, mesh.F_smooth)
        err_smooth = centerline_rms_cs(mesh, N_side)

        # improvement is the primary assertion (ray-starved sparse case
        # has substantial raw error by design); loose absolute ceiling
        @test err_smooth < err_raw
        @test err_smooth < 0.05
        # per-bin energy conservation from the Woodbury solve
        @test maximum(abs.(mesh.energy_error)) < 1e-6
    end
end