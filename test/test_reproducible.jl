println("\n" * "-"^60)
println("Testing reproducibility")
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
    return RayTracingDomain2D([face], [(NDIM_REPRO, NDIM_REPRO)], verbose = false)
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
    # default sampler (Sobol): one integer seed selects the realisation
    a = build_2d_repro(); a(RAYS_REPRO; method = :exchange, verbose = false)
    b = build_2d_repro(); b(RAYS_REPRO; method = :exchange, verbose = false)
    @test a.F_raw == b.F_raw

    c = build_2d_repro(); c(RAYS_REPRO; method = :exchange, seeds = 2, verbose = false)
    @test c.F_raw != a.F_raw

    d = build_2d_repro(); d(RAYS_REPRO; method = :exchange, nthreads = 1, verbose = false)
    @test d.F_raw == a.F_raw                                   # independent of the number of threads

    # pseudorandom sampler: one seed per thread
    ar = build_2d_repro(); ar(RAYS_REPRO; method = :exchange, sampler = :random, verbose = false)
    br = build_2d_repro(); br(RAYS_REPRO; method = :exchange, sampler = :random, verbose = false)
    @test ar.F_raw == br.F_raw

    cr = build_2d_repro()
    cr(RAYS_REPRO; method = :exchange, sampler = :random, seeds = shifted_seeds(), verbose = false)
    @test cr.F_raw != ar.F_raw
end

# ---- 3D view factors ---------------------------------------------------------

@testset "3D view factor determinism" begin
    a = build_3d_vf(); a(; parallel = true, verbose = false)
    b = build_3d_vf(); b(; parallel = true, verbose = false)
    @test a.F_raw == b.F_raw
end

# ---- 3D ray tracing ----------------------------------------------------------

@testset "3D ray tracing reproducibility" begin
    # default sampler (Sobol): one integer seed selects the realisation
    a = build_3d_mc(); a(RAYS_REPRO; verbose = false)
    b = build_3d_mc(); b(RAYS_REPRO; verbose = false)
    @test a.F_raw == b.F_raw

    c = build_3d_mc(); c(RAYS_REPRO; seeds = 2, verbose = false)
    @test c.F_raw != a.F_raw

    e = build_3d_mc(); e(RAYS_REPRO; nthreads = 1, verbose = false)
    @test e.F_raw == a.F_raw                                   # independent of the number of threads

    # pseudorandom sampler: one seed per thread
    ar = build_3d_mc(); ar(RAYS_REPRO; sampler = :random, verbose = false)
    br = build_3d_mc(); br(RAYS_REPRO; sampler = :random, verbose = false)
    @test ar.F_raw == br.F_raw
    cr = build_3d_mc(); cr(RAYS_REPRO; sampler = :random, seeds = shifted_seeds(), verbose = false)
    @test cr.F_raw != ar.F_raw

    # RNG state must not carry over between calls on the same domain
    d = build_3d_mc()
    d(RAYS_REPRO; verbose = false); F1 = copy(d.F_raw)
    d(RAYS_REPRO; verbose = false); F2 = copy(d.F_raw)
    @test F1 == F2
end

# ---- integer seeds -----------------------------------------------------------

@testset "integer seeds select disjoint blocks" begin
    nt = Threads.nthreads()
    # pseudorandom sampler: an integer seed is a block of nthreads per-thread seeds
    default = build_2d_repro(); default(RAYS_REPRO; method = :exchange, sampler = :random, verbose = false)
    one = build_2d_repro(); one(RAYS_REPRO; method = :exchange, sampler = :random, seeds = 1, verbose = false)
    @test one.F_raw == default.F_raw                           # seeds = 1 is the default 1:nthreads
    two = build_2d_repro(); two(RAYS_REPRO; method = :exchange, sampler = :random, seeds = 2, verbose = false)
    blk = build_2d_repro(); blk(RAYS_REPRO; method = :exchange, sampler = :random, seeds = (nt + 1):(2nt), verbose = false)
    @test two.F_raw == blk.F_raw                               # seeds = 2 is the next block
    @test two.F_raw != one.F_raw
    @test_throws ErrorException build_2d_repro()(RAYS_REPRO; method = :exchange, sampler = :random, seeds = 0, verbose = false)

    # default sampler (Sobol): an integer seed is the realisation; per-thread seeds are rejected
    sdef = build_2d_repro(); sdef(RAYS_REPRO; method = :exchange, verbose = false)
    sone = build_2d_repro(); sone(RAYS_REPRO; method = :exchange, seeds = 1, verbose = false)
    @test sone.F_raw == sdef.F_raw                             # seeds = 1 is the default
    @test_throws ErrorException build_2d_repro()(RAYS_REPRO; method = :exchange, seeds = 0, verbose = false)
    @test_throws ErrorException build_2d_repro()(RAYS_REPRO; method = :exchange, seeds = 1:nt, verbose = false)
    
    # 3D, default sampler (Sobol): seeds = 1 is the default realisation
    a = build_3d_mc(); a(RAYS_REPRO; verbose = false)
    b = build_3d_mc(); b(RAYS_REPRO; seeds = 1, verbose = false)
    c = build_3d_mc(); c(RAYS_REPRO; seeds = 2, verbose = false)
    @test a.F_raw == b.F_raw
    @test c.F_raw != a.F_raw
    @test_throws ErrorException build_3d_mc()(RAYS_REPRO; seeds = 1:nt, verbose = false)

    # 3D, pseudorandom sampler: integer seeds select blocks of per-thread seeds
    ar = build_3d_mc(); ar(RAYS_REPRO; sampler = :random, verbose = false)
    b1 = build_3d_mc(); b1(RAYS_REPRO; sampler = :random, seeds = 1, verbose = false)
    b2 = build_3d_mc(); b2(RAYS_REPRO; sampler = :random, seeds = (nt + 1):(2nt), verbose = false)
    c2 = build_3d_mc(); c2(RAYS_REPRO; sampler = :random, seeds = 2, verbose = false)
    @test ar.F_raw == b1.F_raw
    @test c2.F_raw == b2.F_raw
end

@testset "the domain constructor leaves its input faces untouched" begin
    verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    face = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), 1, 1.0, 0.0)
    face.T_in_w  = [1000.0, 0.0, 0.0, 0.0]              # boundary values
    face.epsilon = [1.0, 1.0, 1.0, 1.0]                 # black walls
    face.T_in_g  = -1.0                                 # gas temperature unknown
    face.q_in_g  = 0.0                                  # radiative equilibrium
    a = RayTracingDomain2D([face], [(3, 3)], verbose = false)
    @test isempty(face.subVolumes)                      # the input was not meshed
    @test a.coarse_mesh[1] !== face                     # the domain owns a copy
    b = RayTracingDomain2D([face], [(4, 4)], verbose = false)   # the same face, reused
    @test length(a.fine_mesh[1]) == 9                   # the first domain is unaffected
    @test length(b.fine_mesh[1]) == 16                  # the second one meshed correctly
end

println("✓ Reproducibility tests complete")