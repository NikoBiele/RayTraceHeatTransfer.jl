println("\n" * "-"^60)
println("Testing Sobol quasi random sampler")
println("-"^60)

using Test
using RayTraceHeatTransfer
using GeometryBasics, StaticArrays, SparseArrays
using Random, StatsBase, LinearAlgebra

const RT   = RayTraceHeatTransfer                       # shorthand for unexported internals
const RAYS = 32 * 4096                                  # 4096 rays per emitter on the 4x4 cavity

# ------------------------------------------------------------------ domains ---

# grey unit cavity, black walls, 4x4 cells: 16 wall elements + 16 volumes = 32 emitters
function build_cavity()
    verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    face = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), 1, 1.0, 0.0)
    face.T_in_w  = [1000.0, 0.0, 0.0, 0.0]              # boundary values, unused by the trace
    face.epsilon = [1.0, 1.0, 1.0, 1.0]                 # black walls
    face.T_in_g  = -1.0                                 # gas temperature unknown
    face.q_in_g  = 0.0                                  # radiative equilibrium
    return RayTracingDomain2D([face], [(4, 4)], verbose = false)
end

# circular domain of triangular coarse faces, isothermal rim (from test_triangle_mesh.jl)
function build_triangle_disc(; N_seg = 8, R = 1.0, T_hot = 1000.0, kappa = 1.0)
    rim(j) = Point2(R * cos(2π * (j - 1) / N_seg), R * sin(2π * (j - 1) / N_seg))   # rim vertex j
    faces     = PolyVolume2D{Float64}[]                 # one triangle per segment
    divisions = Tuple{Int,Int}[]                        # subdivision of every triangle
    for j in 1:N_seg
        verts = SVector(Point2(0.0, 0.0), rim(j), rim(j + 1))        # centre, rim, rim
        face = PolyVolume2D{Float64}(verts, SVector(false, true, false), 1, kappa, 0.0)  # spoke, rim, spoke
        face.T_in_w  = [0.0, T_hot, 0.0]                # only the rim wall is real
        face.epsilon = [1.0, 1.0, 1.0]                  # black rim
        face.T_in_g  = -1.0                             # gas temperature unknown
        face.q_in_g  = 0.0                              # radiative equilibrium
        push!(faces, face)
        push!(divisions, (2, 2))
    end
    return RayTracingDomain2D(faces, divisions, verbose = false)
end

# two square zones side by side, two spectral bins: bin 1 has the same kappa in both
# zones (uniform), bin 2 differs between the zones (nonuniform)
function build_two_zone()
    function zone(x0, solid, kappa)                     # unit square starting at x = x0
        n = length(kappa)                               # number of spectral bins
        verts = SVector(Point2(x0, 0.0), Point2(x0 + 1.0, 0.0), Point2(x0 + 1.0, 1.0), Point2(x0, 1.0))
        face = PolyVolume2D{Float64}(verts, solid, n, 1.0, 0.0)
        face.kappa_g   = copy(kappa)                    # absorption coefficient per bin
        face.sigma_s_g = zeros(n)                       # no scattering
        face.epsilon   = [fill(1.0, n) for _ in 1:4]    # black walls, per wall and bin
        face.T_in_w    = fill(1000.0, 4)                # wall temperatures, unused by the trace
        face.q_in_w    = zeros(4)                       # wall fluxes, unused by the trace
        face.T_in_g    = -1.0                           # gas temperature unknown
        face.q_in_g    = 0.0                            # radiative equilibrium
        return face
    end
    left  = zone(0.0, SVector(true, false, true, true), [0.1, 2.0])   # right wall (wall 2) open
    right = zone(1.0, SVector(true, true, true, false), [0.1, 5.0])   # left wall (wall 4) open
    return RayTracingDomain2D([left, right], [(3, 3), (3, 3)], verbose = false)
end

# ------------------------------------------------------------------ helpers ---

# trace the cavity with the given method and keywords; dense copy of F_raw
function trace_cavity(method::Symbol; rays = RAYS, kwargs...)
    mesh = build_cavity()                               # fresh domain for every trace
    mesh(rays; method = method, verbose = false, kwargs...)
    return Matrix(mesh.F_raw)
end

# trace the two-zone domain; one dense matrix per spectral bin
function trace_two_zone(method::Symbol; rays = 2_000_000, kwargs...)
    mesh = build_two_zone()                             # fresh domain for every trace
    mesh(rays; method = method, verbose = false, kwargs...)
    return [Matrix(F) for F in mesh.F_raw]
end

maxdiff(A, B) = maximum(abs.(A .- B))                   # largest entry-wise difference
entry_noise(Fs) = sqrt(mean(var(cat(Fs...; dims = 3); dims = 3)))   # rms spread of an entry across realisations

# best of n timings of one cavity trace (after a warm-up, so compilation is not timed)
function best_time(method::Symbol, sampler::Symbol; rays = 32 * 65536, n = 3)
    trace_cavity(method; rays = 32 * 256, sampler = sampler)        # warm-up
    return minimum(@elapsed(trace_cavity(method; rays = rays, sampler = sampler)) for _ in 1:n)
end

# ------------------------------------------------------------------- checks ---

@testset "Sobol sampler integration" begin

    @testset "sampler: points, shifts, contract" begin
        # first n rays of one emitter as an n x d matrix
        function first_points(n, d; seed = 1, bin = 1, emitter = 7, T = Float64)
            source = RT.SobolEmitterRNG(d, seed, bin, emitter)      # one emitter's source
            pts = zeros(T, n, d)                                    # one row per ray (not named X: see below)
            for i in 1:n
                RT.nextRay!(source)                                 # announce the ray
                for j in 1:d
                    pts[i, j] = rand(source, T)                     # coordinates in order
                end
            end
            return pts                                              # an inner function shares the testset's locals by name
        end
        X = first_points(1024, 5)
        # every coordinate puts exactly one of the first 1024 points in each of 1024 intervals
        @test all(sort(floor.(Int, X[:, j] .* 1024)) == collect(0:1023) for j in 1:5)
        @test first_points(1024, 5) == X                            # reproducible
        @test first_points(1024, 5; emitter = 8) != X               # shift depends on the emitter
        @test first_points(1024, 5; seed = 2) != X                  # ... on the seed
        @test first_points(1024, 5; bin = 2) != X                   # ... and on the spectral bin
        X32 = first_points(1024, 5; T = Float32)
        @test all(0 .<= X32 .< 1)                                   # Float32 draws stay inside [0, 1)
        # the leading coordinates do not depend on the total dimension (beta = 0 draws one fewer)
        @test first_points(256, 4)[:, 1:3] == first_points(256, 3)
        # numbers per ray: wall 3, triangle 4, quad 5, plus one when the depth is sampled
        @test RT.sobolDims2D(true, 4, false) == 3 && RT.sobolDims2D(true, 4, true) == 4
        @test RT.sobolDims2D(false, 3, false) == 4 && RT.sobolDims2D(false, 3, true) == 5
        @test RT.sobolDims2D(false, 4, false) == 5 && RT.sobolDims2D(false, 4, true) == 6
        # drawing before nextRay!, or more numbers than the ray owns, is an error
        q = RT.SobolEmitterRNG(3, 1, 1, 7)
        @test_throws "draw order contract" rand(q, Float64)
        RT.nextRay!(q); rand(q, Float64); rand(q, Float64); rand(q, Float64)
        @test_throws "draw order contract" rand(q, Float64)
    end

    @testset ":pathlength with the Sobol default" begin
        F = trace_cavity(:pathlength)                               # default sampler, seed 1
        @test maximum(abs.(sum(F, dims = 2) .- 1)) < 1e-12          # rows conserve energy
        @test maxdiff(trace_cavity(:pathlength), F) == 0.0          # reproducible
        @test maxdiff(trace_cavity(:pathlength; nthreads = 1), F) < 1e-12              # independent of threads
        @test maxdiff(trace_cavity(:pathlength; chunk_rays = 32 * 1000), F) < 1e-12    # 5 chunks, same rays
        @test maxdiff(trace_cavity(:pathlength; chunk_rays = 32 * 333), F) < 1e-12     # 13 uneven chunks
        @test maxdiff(trace_cavity(:pathlength; sampler = :sobol), F) == 0.0           # explicit = default
        @test maxdiff(trace_cavity(:pathlength; seeds = 2), F) > 1e-4                  # another realisation
        @test maxdiff(trace_cavity(:pathlength; sampler = :random), F) > 1e-4          # pseudorandom still available
    end

    @testset ":exchange with the Sobol default" begin
        F = trace_cavity(:exchange)                                 # default sampler, seed 1
        @test maximum(abs.(sum(F, dims = 2) .- 1)) < 1e-12          # rows conserve energy
        @test maxdiff(trace_cavity(:exchange), F) == 0.0            # reproducible
        @test maxdiff(trace_cavity(:exchange; nthreads = 1), F) == 0.0                 # integer counts: exactly equal
        @test maxdiff(trace_cavity(:exchange; sampler = :sobol), F) == 0.0             # explicit = default
        @test maxdiff(trace_cavity(:exchange; seeds = 2), F) > 1e-4                    # another realisation
        @test maxdiff(trace_cavity(:exchange; sampler = :random), F) > 1e-4            # pseudorandom still available
        # pseudorandom sampling keeps its old seed semantics
        Fr = trace_cavity(:exchange; sampler = :random, seeds = 3)
        @test maxdiff(trace_cavity(:exchange; sampler = :random, seeds = 3), Fr) == 0.0
        @test maxdiff(trace_cavity(:exchange; sampler = :random, rngs = Xoshiro(1), seeds = 3), Fr) == 0.0
    end

    @testset "keyword errors" begin
        for method in (:exchange, :pathlength)
            @test_throws "Pass sampler = :random to use your own generators" trace_cavity(method; rngs = Xoshiro(1))
            @test_throws "takes a single integer seed" trace_cavity(method; seeds = 1:Threads.nthreads())
            @test_throws "must be ≥ 1" trace_cavity(method; seeds = 0)
            @test_throws "Unknown sampler" trace_cavity(method; sampler = :halton)
        end
        @test_throws "not available for method = :direct" trace_cavity(:direct; sampler = :sobol)
    end

    @testset "unbiased, and quieter than pseudorandom" begin
        n = 8                                                       # realisations per sampler
        for (method, min_gain) in ((:exchange, 1.2), (:pathlength, 2.0))
            Fs = [trace_cavity(method; seeds = s) for s in 1:n]                        # Sobol realisations
            Fr = [trace_cavity(method; seeds = s, sampler = :random) for s in 1:n]     # pseudorandom realisations
            n_s, n_r = entry_noise(Fs), entry_noise(Fr)             # noise of a single realisation
            bias = sqrt(mean(abs2, mean(Fs) .- mean(Fr)))           # rms difference of the two means
            expected = sqrt(n_s^2 + n_r^2) / sqrt(n)                # what pure noise would give
            @test bias < 2 * expected                               # the two samplers agree on the mean
            @test n_r / n_s > min_gain                              # Sobol is clearly quieter
        end
    end

    @testset "triangular emitters: isothermal disc is exact" begin
        mesh = build_triangle_disc()                                # isothermal rim at 1000 K
        mesh(2_000_000; method = :exchange, verbose = false)        # Sobol default; triangles draw 4 + 1 numbers
        smooth!(mesh, verbose = false)
        solveEquilibrium!(mesh, mesh.F_smooth, verbose = false)
        T_g = [ff.T_g for fine in mesh.fine_mesh for ff in fine]    # all gas temperatures
        @test maximum(abs.(T_g .- 1000.0)) < 1e-3                   # gas must sit at the rim temperature
    end

    @testset "variable spectral :exchange (uniform and nonuniform bins)" begin
        probe = build_two_zone()
        @test probe.uniform_across_bin[1] > 0                       # bin 1 takes the uniform-group path
        @test probe.uniform_across_bin[2] < 0                       # bin 2 takes the nonuniform path
        F = trace_two_zone(:exchange)                               # Sobol default
        @test length(F) == 2
        @test all(maximum(abs.(sum(F[k], dims = 2) .- 1)) < 1e-12 for k in 1:2)       # rows conserve energy
        @test F[1] != F[2]                                          # the bins really differ
        F_again = trace_two_zone(:exchange)
        @test all(maxdiff(F_again[k], F[k]) == 0.0 for k in 1:2)    # reproducible
        F_one = trace_two_zone(:exchange; nthreads = 1)
        @test all(maxdiff(F_one[k], F[k]) == 0.0 for k in 1:2)      # independent of threads
        F_rand = trace_two_zone(:exchange; sampler = :random)
        @test all(norm(F[k] - F_rand[k]) / norm(F_rand[k]) < 5e-2 for k in 1:2)       # agrees with pseudorandom
        F_path = trace_two_zone(:pathlength)
        @test all(norm(F[k] - F_path[k]) / norm(F_path[k]) < 3e-2 for k in 1:2)       # agrees with :pathlength
    end

    @testset "no slowdown from the sampler" begin
        for method in (:exchange, :pathlength)
            t_s = best_time(method, :sobol)                         # Sobol default
            t_r = best_time(method, :random)                        # pseudorandom
            @test t_s < 1.5 * t_r                                   # a boxed generator would show up here
        end
    end

end

const CUBE_POINTS = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0; 0.0 1.0 1.0;
                     1.0 0.0 0.0; 1.0 0.0 1.0; 1.0 1.0 0.0; 1.0 1.0 1.0]      # unit cube corners
const CUBE_FACES  = [1 2 4 3; 5 6 8 7; 1 5 7 3; 2 6 8 4; 3 4 8 7; 1 2 6 5]    # its six faces
const CUBE_TIN    = [1000.0, 0.0, -1.0, -1.0, -1.0, -1.0]                     # boundary values, unused by the trace
const CUBE_QIN    = [-1.0, -1.0, 0.0, 0.0, 0.0, 0.0]                          # boundary values, unused by the trace
const NDIM        = 3                                                         # 3x3 elements per face: 54 elements

# trace the cube with the given keywords; dense copy of F_raw
function trace_cube(; rays_per_emitter = 4096, kwargs...)
    d = RayTracingDomain3D_surfaces(CUBE_POINTS, CUBE_FACES, NDIM, CUBE_QIN, CUBE_TIN, ones(6))
    d(rays_per_emitter * 6 * NDIM^2; verbose = false, kwargs...)
    return Matrix(d.F_raw)
end

# face-to-face view factors from the element matrix (elements ordered face by face, equal areas)
function face_factors(F)
    n = NDIM^2                                           # elements per face
    FF = zeros(6, 6)
    for I in 1:6, J in 1:6
        rows = (I - 1) * n + 1 : I * n                   # elements of face I
        cols = (J - 1) * n + 1 : J * n                   # elements of face J
        FF[I, J] = mean(sum(F[rows, cols], dims = 2))    # mean over emitters of the fraction reaching face J
    end
    return FF
end

const FF_EXACT = [I == J ? 0.0 : 0.2 for I in 1:6, J in 1:6]                  # cube: 1/5 to each other face

face_error(Fs) = sqrt(mean(mean(abs2, face_factors(F) .- FF_EXACT) for F in Fs))   # rms error of the 36 face factors

# best of n timings of one trace (after a warm-up, so compilation is not timed)
function best_time(sampler::Symbol; rays_per_emitter = 65536, n = 3)
    trace_cube(rays_per_emitter = 256, sampler = sampler)                     # warm-up
    return minimum(@elapsed(trace_cube(rays_per_emitter = rays_per_emitter, sampler = sampler)) for _ in 1:n)
end

@testset "Sobol sampler in the 3D surface tracer" begin

    @testset "Sobol default" begin
        F = trace_cube()                                                      # default sampler, seed 1
        @test maximum(abs.(sum(F, dims = 2) .- 1)) < 1e-12                    # closed cube: nothing escapes
        @test maxdiff(trace_cube(), F) == 0.0                                 # reproducible
        @test maxdiff(trace_cube(nthreads = 1), F) == 0.0                     # independent of threads
        @test maxdiff(trace_cube(sampler = :sobol), F) == 0.0                 # explicit = default
        @test maxdiff(trace_cube(seeds = 2), F) > 1e-4                        # another realisation
        @test maxdiff(trace_cube(sampler = :random), F) > 1e-4                # pseudorandom still available
    end

    @testset "pseudorandom sampler keeps its seed contract" begin
        Fr = trace_cube(sampler = :random, seeds = 3)
        @test maxdiff(trace_cube(sampler = :random, seeds = 3), Fr) == 0.0
        @test maxdiff(trace_cube(sampler = :random, seeds = 3, rngs = Xoshiro(1)), Fr) == 0.0
        @test maxdiff(trace_cube(sampler = :random, nthreads = 1), trace_cube(sampler = :random, nthreads = 1)) == 0.0
    end

    @testset "keyword errors" begin
        @test_throws "Pass sampler = :random to use your own generators" trace_cube(rngs = Xoshiro(1))
        @test_throws "takes a single integer seed" trace_cube(seeds = 1:Threads.nthreads())
        @test_throws "must be ≥ 1" trace_cube(seeds = 0)
        @test_throws "Unknown sampler" trace_cube(sampler = :halton)
    end

    @testset "unbiased, more accurate and quieter than pseudorandom" begin
        n = 8                                                                 # realisations per sampler
        Fs = [trace_cube(seeds = s) for s in 1:n]                             # Sobol realisations
        Fr = [trace_cube(seeds = s, sampler = :random) for s in 1:n]          # pseudorandom realisations
        n_s, n_r = entry_noise(Fs), entry_noise(Fr)                           # noise of a single realisation
        e_s, e_r = face_error(Fs), face_error(Fr)                             # error against the exact face factors
        bias = sqrt(mean(abs2, mean(Fs) .- mean(Fr)))                         # rms difference of the two means
        expected = sqrt(n_s^2 + n_r^2) / sqrt(n)                              # what pure noise would give
        @test bias < 2 * expected                                             # the two samplers agree on the mean
        @test n_r / n_s > 1.2                                                 # quieter entry by entry
        @test e_r / e_s > 1.5                                                 # and more accurate against the exact answer
    end

    @testset "no slowdown from the sampler" begin
        t_s = best_time(:sobol)                                               # Sobol default
        t_r = best_time(:random)                                              # pseudorandom
        @test t_s < 1.5 * t_r
    end

end

println("✓ Sobol quasi random sampling tests complete")