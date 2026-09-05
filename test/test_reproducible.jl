println("\n" * "-"^60)
println("Testing Reproducibility")
println("-"^60)

using RayTraceHeatTransfer
using Test
using StaticArrays
using GeometryBasics
using SparseArrays

const NDIM_REPRO = 3
const RAYS_REPRO = 50_000

# ---- geometry builders -------------------------------------------------------

function build_2d_repro()
    vertices   = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0),
                         Point2(1.0, 1.0), Point2(0.0, 1.0))
    solidWalls = SVector(true, true, true, true)
    face = PolyVolume2D{Float64}(vertices, solidWalls, 1, 1.0, 0.0)
    face.T_in_w  = [1000.0, 0.0, 0.0, 0.0]
    face.epsilon = [1.0, 1.0, 1.0, 1.0]
    face.T_in_g  = -1.0
    face.q_in_g  = 0.0
    return RayTracingDomain2D([face], [(NDIM_REPRO, NDIM_REPRO)])
end

const CUBE_POINTS = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0; 0.0 1.0 1.0;
                     1.0 0.0 0.0; 1.0 0.0 1.0; 1.0 1.0 0.0; 1.0 1.0 1.0]
const CUBE_FACES  = [1 2 4 3; 5 6 8 7; 1 5 7 3; 2 6 8 4; 3 4 8 7; 1 2 6 5]
const CUBE_TIN    = [1000.0, 0.0, -1.0, -1.0, -1.0, -1.0]
const CUBE_QIN    = [-1.0, -1.0, 0.0, 0.0, 0.0, 0.0]
const CUBE_EPS    = ones(6)

build_3d_mc() = RayTracingDomain3D_surfaces(CUBE_POINTS, CUBE_FACES, NDIM_REPRO,
                                            CUBE_QIN, CUBE_TIN, CUBE_EPS)
build_3d_vf() = ViewFactorDomain3D(CUBE_POINTS, CUBE_FACES, NDIM_REPRO,
                                   CUBE_QIN, CUBE_TIN, CUBE_EPS)

shifted_seeds() = (1:Threads.nthreads()) .+ 1000

# ---- 2D exchange -------------------------------------------------------------

@testset "2D exchange reproducibility" begin
    a = build_2d_repro(); a(RAYS_REPRO; method = :exchange, verbose = false)
    b = build_2d_repro(); b(RAYS_REPRO; method = :exchange, verbose = false)
    @test a.F_raw == b.F_raw

    c = build_2d_repro()
    c(RAYS_REPRO; method = :exchange, seeds = shifted_seeds(), verbose = false)
    @test c.F_raw != a.F_raw
end

# ---- 3D view factors ---------------------------------------------------------

@testset "3D view factor determinism" begin
    a = build_3d_vf(); a(; parallel = true, verbose = false)
    b = build_3d_vf(); b(; parallel = true, verbose = false)
    @test a.F_raw == b.F_raw
end

# ---- 3D ray tracing ----------------------------------------------------------

@testset "3D ray tracing reproducibility" begin
    a = build_3d_mc(); a(RAYS_REPRO; verbose = false)
    b = build_3d_mc(); b(RAYS_REPRO; verbose = false)
    @test a.F_raw == b.F_raw

    c = build_3d_mc(); c(RAYS_REPRO; seeds = shifted_seeds(), verbose = false)
    @test c.F_raw != a.F_raw

    # RNG state must not carry over between calls on the same domain
    d = build_3d_mc()
    d(RAYS_REPRO; verbose = false); F1 = copy(d.F_raw)
    d(RAYS_REPRO; verbose = false); F2 = copy(d.F_raw)
    @test F1 == F2
end