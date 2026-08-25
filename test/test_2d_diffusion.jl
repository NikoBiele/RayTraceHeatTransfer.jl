println("\n" * "-"^60)
println("Testing 2D grey diffusion approximation")
println("-"^60)

using RayTraceHeatTransfer
using Test
using LinearAlgebra
using StatsBase
using StaticArrays
using GeometryBasics
using SparseArrays
using ConvolutionInterpolations

const T_HOT_DIFF   = 1000.0
const BETA_DIFF    = 25.0
const ASPECT_DIFF  = 1000.0

# ---- analytical diffusion reference -----------------------------------------
function diffusion_S(z, beta, D, eps_w1, eps_w2, E_bw1, E_bw2)
    q_z  = (E_bw1 - E_bw2) / (3 * beta * D / 4 + 1 / eps_w1 + 1 / eps_w2 - 1)
    E_b1 = E_bw1 + q_z * (1 / 2 - 1 / eps_w1)
    return E_b1 - (3 * beta * z / 4) * q_z
end

centerline_taus(N_side) = collect(range(1 / (2N_side), 1 - 1 / (2N_side), length = N_side))

# ---- problem builders --------------------------------------------------------
function build_diffusion_grey(N_side::Int)
    vertices = SVector(
        Point2(0.0, 0.0), Point2(ASPECT_DIFF, 0.0),
        Point2(ASPECT_DIFF, 1.0), Point2(0.0, 1.0))
    solidWalls = SVector(true, true, true, true)
    face = PolyVolume2D{Float64}(vertices, solidWalls, 1, BETA_DIFF, 0.0)
    face.T_in_w  = [T_HOT_DIFF, 0.0, 0.0, 0.0]
    face.epsilon = [1.0, 1.0, 1.0, 1.0]
    face.T_in_g  = -1.0
    face.q_in_g  = 0.0
    return RayTracingDomain2D([face], [(N_side, N_side)])
end

# ---- extraction + error ------------------------------------------------------
function centerline_rms_S(mesh, N_side)
    all_T = [ff.T_g for ff in mesh.fine_mesh[1]]
    T_center = reshape(all_T, N_side, N_side)[(N_side - 1) ÷ 2 + 1, :]
    S_comp = (T_center ./ T_HOT_DIFF) .^ 4
    S_ref  = diffusion_S.(centerline_taus(N_side), BETA_DIFF, 1.0, 1.0, 1.0, 1.0, 0.0)
    return sqrt(mean((S_comp .- S_ref) .^ 2))
end

# =============================================================================
@testset "Diffusion-limit validation (dense/sparse x grey)" begin

    # tolerances: statistical (MC ray trace) + diffusion-limit model error
    RMS_TOL   = 0.02     # absolute RMS in S(tau) for the smoothed solve
    RATIO_TOL = 0.5      # smoothed error must be well below raw error

    (label, N_side, expect_sparse) = ("sparse", 31, true)
    N_rays = (4 * N_side + N_side^2) * 1_000

    @testset "$label grey" begin
        mesh = build_diffusion_grey(N_side)
        mesh(N_rays; method = :exchange)
        smooth!(mesh)

        @test (mesh.F_smooth isa SparseMatrixCSC) == expect_sparse
        @test count(<(0.0), mesh.F_smooth) == 0

        # raw solve
        solveEquilibrium!(mesh, mesh.F_raw)
        err_raw = centerline_rms_S(mesh, N_side)

        # smoothed solve
        solveEquilibrium!(mesh, mesh.F_smooth)
        err_ap = centerline_rms_S(mesh, N_side)

        @test err_ap < RMS_TOL
        @test err_ap < RATIO_TOL * err_raw
        @test abs(sum(mesh.energy_error)) < 1e-6
    end
end