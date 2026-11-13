println("\n" * "-"^60)
println("Testing ray allocation over the emitters")
println("-"^60)

using RayTraceHeatTransfer
using Test
using StaticArrays
using GeometryBasics
using SparseArrays
using StatsBase
using RayTraceHeatTransfer: emitter_ray_counts, partition_emitters, get_w

const NDIM_ALLOC = 5
const RAYS_ALLOC = 9_000

# 2 x 1 enclosure with absorption coefficient beta and black walls; the wall
# elements on the long sides are twice as long as those on the short sides
function build_2d_alloc(beta::Real; surfaces_only::Bool = false)
    vertices   = SVector(Point2(0.0, 0.0), Point2(2.0, 0.0),
                         Point2(2.0, 1.0), Point2(0.0, 1.0))
    solidWalls = SVector(true, true, true, true)
    face = PolyVolume2D{Float64}(vertices, solidWalls, 1, beta, 0.0)
    face.T_in_w  = [1000.0, 0.0, 0.0, 0.0]
    face.epsilon = [1.0, 1.0, 1.0, 1.0]
    face.T_in_g  = -1.0
    face.q_in_g  = 0.0
    return RayTracingDomain2D([face], [(NDIM_ALLOC, NDIM_ALLOC)];
                              verbose = false, surfaces_only = surfaces_only)
end

# largest relative spread of counts per unit weight over a set of emitters
function spread_per_weight(counts, w, idx)
    ratio = counts[idx] ./ w[idx]
    return maximum(abs.(ratio .- mean(ratio))) / mean(ratio)
end

@testset "uniform allocation" begin
    println("  uniform allocation")
    mesh = build_2d_alloc(1.0)
    N = length(get_w(mesh))
    counts = emitter_ray_counts(mesh, RAYS_ALLOC, N, :uniform, 1)
    @test all(counts .== div(RAYS_ALLOC, N))
    # equal counts keep the plain split by emitter number
    @test partition_emitters(fill(7, 10), 3) == [1:4, 5:7, 8:10]
end

@testset "weighted allocation, optically thin" begin
    println("  weighted allocation, optically thin")
    mesh = build_2d_alloc(0.1)
    w = get_w(mesh)
    N = length(w)
    n_surf = length(mesh.surface_mapping)
    counts = emitter_ray_counts(mesh, RAYS_ALLOC, N, :weighted, 1)
    @test abs(sum(counts) - RAYS_ALLOC) <= N                    # the budget is kept up to rounding
    # the w split gives the gas less than its equal-count share, so plain w applies
    share_w = sum(w[(n_surf + 1):end]) / sum(w)
    @test share_w < (N - n_surf) / N
    @test abs(sum(counts[(n_surf + 1):end]) / sum(counts) - share_w) < 0.01
    @test spread_per_weight(counts, w, 1:n_surf) < 0.01         # walls in proportion to their length
end

@testset "weighted allocation, optically thick: gas share capped" begin
    println("  weighted allocation, optically thick")
    mesh = build_2d_alloc(20.0)
    w = get_w(mesh)
    N = length(w)
    n_surf = length(mesh.surface_mapping)
    counts = emitter_ray_counts(mesh, RAYS_ALLOC, N, :weighted, 1)
    @test abs(sum(counts) - RAYS_ALLOC) <= N
    # the w split would give the gas more than its equal-count share, so the share is capped
    @test sum(w[(n_surf + 1):end]) / sum(w) > (N - n_surf) / N
    @test abs(sum(counts[(n_surf + 1):end]) / sum(counts) - (N - n_surf) / N) < 0.01
    @test spread_per_weight(counts, w, 1:n_surf) < 0.01         # walls still in proportion to their length
    @test spread_per_weight(counts, w, (n_surf + 1):N) < 0.01   # gas in proportion to its weight
end

@testset "surfaces-only domain: one ray per gas cell" begin
    println("  surfaces-only domain")
    mesh = build_2d_alloc(5.0; surfaces_only = true)
    w = get_w(mesh)
    N = length(w)
    n_surf = length(mesh.surface_mapping)
    counts = emitter_ray_counts(mesh, RAYS_ALLOC, N, :weighted, 1)
    @test all(counts[(n_surf + 1):end] .== 1)
    @test spread_per_weight(counts, w, 1:n_surf) < 0.01
end

@testset "thread partition balances rays" begin
    println("  thread partition")
    counts = vcat([400], fill(10, 60), [250, 250])
    n = length(counts)
    nt = 4
    ranges = partition_emitters(counts, nt)
    @test length(ranges) == nt
    @test reduce(vcat, collect.(ranges)) == collect(1:n)        # contiguous and complete
    loads = [sum(counts[r]; init = 0) for r in ranges]
    @test maximum(loads) <= sum(counts) / nt + maximum(counts)
end

@testset "weighted trace" begin
    println("  weighted trace")
    a = build_2d_alloc(20.0); a(RAYS_ALLOC; method = :exchange, verbose = false)
    @test maximum(abs.(vec(sum(a.F_raw; dims = 2)) .- 1)) < 1e-12
    b = build_2d_alloc(20.0); b(RAYS_ALLOC; method = :exchange, verbose = false)
    @test a.F_raw == b.F_raw                                     # reproducible
    c = build_2d_alloc(20.0); c(RAYS_ALLOC; method = :exchange, nthreads = 1, verbose = false)
    @test c.F_raw == a.F_raw                                     # independent of the number of threads
    d = build_2d_alloc(20.0); d(RAYS_ALLOC; method = :exchange, ray_allocation = :uniform, verbose = false)
    @test d.F_raw != a.F_raw                                     # the allocation changes the realisation
end

@testset "weighted equals uniform when all weights are equal" begin
    println("  weighted equals uniform for equal weights")
    # unit square on an n x n mesh: walls 1/n long, gas cells 1/n^2 in area, so
    # 4 beta / n^2 = 1/n, i.e. beta = n/4, gives every emitter the same weight
    function build_square(beta::Real)
        vertices   = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0),
                             Point2(1.0, 1.0), Point2(0.0, 1.0))
        solidWalls = SVector(true, true, true, true)
        face = PolyVolume2D{Float64}(vertices, solidWalls, 1, beta, 0.0)
        face.T_in_w  = [1000.0, 0.0, 0.0, 0.0]
        face.epsilon = [1.0, 1.0, 1.0, 1.0]
        face.T_in_g  = -1.0
        face.q_in_g  = 0.0
        return RayTracingDomain2D([face], [(NDIM_ALLOC, NDIM_ALLOC)]; verbose = false)
    end
    beta_equal = NDIM_ALLOC / 4
    sq = build_square(beta_equal)
    w = get_w(sq)
    N = length(w)
    @test RAYS_ALLOC % N == 0                                    # div and round give the same count
    @test maximum(abs.(w ./ w[1] .- 1)) < 1e-12                  # every emitter carries the same weight
    @test emitter_ray_counts(sq, RAYS_ALLOC, N, :weighted, 1) == emitter_ray_counts(sq, RAYS_ALLOC, N, :uniform, 1)
    u = build_square(beta_equal); u(RAYS_ALLOC; method = :exchange, ray_allocation = :uniform, verbose = false)
    v = build_square(beta_equal); v(RAYS_ALLOC; method = :exchange, ray_allocation = :weighted, verbose = false)
    @test u.F_raw == v.F_raw                                     # same seed, same counts: the same trace
end

@testset "allocation keyword checks" begin
    println("  allocation keyword checks")
    mesh = build_2d_alloc(1.0)
    @test_throws ErrorException mesh(RAYS_ALLOC; method = :pathlength, ray_allocation = :weighted, verbose = false)
    @test_throws ErrorException mesh(RAYS_ALLOC; ray_allocation = :proportional, verbose = false)
end