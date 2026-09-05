println("\n" * "-"^60)
println("Testing 3D Surface Ray Tracing")
println("-"^60)

using RayTraceHeatTransfer
using Test
using LinearAlgebra
using SparseArrays
using StatsBase

# ---- geometry ----------------------------------------------------------------

const CUBE_POINTS = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0; 0.0 1.0 1.0;
                     1.0 0.0 0.0; 1.0 0.0 1.0; 1.0 1.0 0.0; 1.0 1.0 1.0]
const CUBE_FACES  = [1 2 4 3; 5 6 8 7; 1 5 7 3; 2 6 8 4; 3 4 8 7; 1 2 6 5]
const CUBE_TIN    = [1000.0, 0.0, -1.0, -1.0, -1.0, -1.0]
const CUBE_QIN    = [-1.0, -1.0, 0.0, 0.0, 0.0, 0.0]

# regular tetrahedron, volume 8/3
const TET_POINTS = [ 1.0  1.0  1.0;  1.0 -1.0 -1.0; -1.0  1.0 -1.0; -1.0 -1.0  1.0]
const TET_FACES  = [1 2 3; 1 2 4; 1 3 4; 2 3 4]
const TET_TIN    = [1000.0, 0.0, -1.0, -1.0]
const TET_QIN    = [-1.0, -1.0, 0.0, 0.0]

# L-shaped duct: 8 cross-section vertices extruded z = 0 -> 1, volume 3
const L_XY = [0.0 0.0; 1.0 0.0; 2.0 0.0; 2.0 1.0; 1.0 1.0; 0.0 1.0; 1.0 2.0; 0.0 2.0]
const L_POINTS = vcat(hcat(L_XY, zeros(8)), hcat(L_XY, ones(8)))
const L_FACES  = [ 1  2  5  6;  2  3  4  5;  6  5  7  8;      # z = 0 caps
                   9 10 13 14; 10 11 12 13; 14 13 15 16;      # z = 1 caps
                   1  2 10  9;  2  3 11 10;  3  4 12 11;      # side walls
                   4  5 13 12;  5  7 15 13;  7  8 16 15;
                   8  6 14 16;  6  1  9 14]
const L_TIN = [1000.0, 0.0, fill(-1.0, 12)...]
const L_QIN = [  -1.0, -1.0, zeros(12)...]

# ---- helpers -----------------------------------------------------------------

collectT(d) = [sf.T_w for f in d.facesMesh for sf in f.subFaces]
collectTin(d) = [sf.T_in_w for f in d.facesMesh for sf in f.subFaces]

function faceRanges(d)
    r = UnitRange{Int}[]
    i = 1
    for f in d.facesMesh
        push!(r, i:(i + length(f.subFaces) - 1))
        i += length(f.subFaces)
    end
    return r
end

faceMeans(d) = [mean(collectT(d)[r]) for r in faceRanges(d)]

function tracedDomain(points, faces, Ndim, q_in, T_in, eps, rays_per_emitter)
    d = RayTracingDomain3D_surfaces(points, faces, Ndim, q_in, T_in, eps)
    d(rays_per_emitter * length(d.facesMesh) * Ndim^2; verbose = false)
    return d
end

function solvedMC(points, faces, Ndim, q_in, T_in, eps, rays_per_emitter)
    d = tracedDomain(points, faces, Ndim, q_in, T_in, eps, rays_per_emitter)
    smooth!(d; verbose = false)
    solveEquilibrium!(d, d.F_smooth; verbose = false)
    return d
end

function solvedVF(points, faces, Ndim, q_in, T_in, eps)
    d = ViewFactorDomain3D(points, faces, Ndim, q_in, T_in, eps)
    d(; parallel = true, verbose = false)
    smooth!(d; verbose = false)
    solveEquilibrium!(d, d.F_smooth; verbose = false)
    return d
end

#############################################################################
### TEST 1: WINDING INDEPENDENCE ############################################
#############################################################################

@testset "Winding independence" begin
    Ndim = 3
    rpe  = 200_000
    eps6 = ones(6)

    ref = solvedVF(CUBE_POINTS, CUBE_FACES, Ndim, CUBE_QIN, CUBE_TIN, eps6)
    ref_means = faceMeans(ref)

    variants = Dict(
        "as given"      => CUBE_FACES,
        "all reversed"  => CUBE_FACES[:, end:-1:1],
        "rotated rows"  => vcat([circshift(CUBE_FACES[i, :], i)' for i in 1:6]...))

    for (label, F) in variants
        @testset "$label" begin
            d = solvedMC(CUBE_POINTS, F, Ndim, CUBE_QIN, CUBE_TIN, eps6, rpe)
            @test abs(d.energy_error) < 1e-10
            @test maximum(abs, faceMeans(d) .- ref_means) < 2.0
        end
    end
end

#############################################################################
### TEST 2: BROKEN GEOMETRY IS REJECTED #####################################
#############################################################################

@testset "Broken geometry rejected at construction" begin
    eps5 = ones(5); eps7 = ones(7)

    holed = CUBE_FACES[2:end, :]                       # one face missing
    @test_throws ErrorException RayTracingDomain3D_surfaces(
        CUBE_POINTS, holed, 3, CUBE_QIN[2:end], CUBE_TIN[2:end], eps5)

    doubled = vcat(CUBE_FACES, CUBE_FACES[1:1, :])     # one face duplicated
    @test_throws ErrorException RayTracingDomain3D_surfaces(
        CUBE_POINTS, doubled, 3, [CUBE_QIN; -1.0], [CUBE_TIN; 1000.0], eps7)
end

#############################################################################
### TEST 3: TRACER INVARIANTS ###############################################
#############################################################################

@testset "Tracer invariants" begin
    cases = [("cube",        CUBE_POINTS, CUBE_FACES, CUBE_QIN, CUBE_TIN, ones(6)),
             ("tetrahedron", TET_POINTS,  TET_FACES,  TET_QIN,  TET_TIN,  ones(4)),
             ("L-duct",      L_POINTS,    L_FACES,    L_QIN,    L_TIN,    ones(14))]

    for (name, pts, fcs, q, T, e) in cases
        @testset "$name" begin
            d = tracedDomain(pts, fcs, 3, q, T, e, 50_000)
            F = d.F_raw

            @test maximum(abs, diag(Matrix(F))) == 0.0        # cannot see itself
            @test maximum(abs, 1 .- vec(sum(F, dims = 2))) < 1e-12
            @test minimum(F) >= 0.0
        end
    end
end

#############################################################################
### TEST 4: AGREEMENT WITH ANALYTICAL VIEW FACTORS ##########################
#############################################################################

@testset "MC vs analytical view factors" begin
    Ndim = 3
    rpe  = 500_000

    cases = [("cube",        CUBE_POINTS, CUBE_FACES, CUBE_QIN, CUBE_TIN, ones(6)),
             ("tetrahedron", TET_POINTS,  TET_FACES,  TET_QIN,  TET_TIN,  ones(4))]

    for (name, pts, fcs, q, T, e) in cases
        @testset "$name" begin
            mc = tracedDomain(pts, fcs, Ndim, q, T, e, rpe)
            vf = ViewFactorDomain3D(pts, fcs, Ndim, q, T, e)
            vf(; parallel = true, verbose = false)

            FMC = Matrix(mc.F_raw)
            FVF = Matrix(vf.F_raw)
            @test size(FMC) == size(FVF)

            # per-entry binomial standard deviation at rpe samples per row
            sd = sqrt.(max.(FVF .* (1 .- FVF), 0.0) ./ rpe)
            z  = maximum(abs.(FMC .- FVF) ./ max.(sd, 1e-12))
            @test z < 6.0
        end
    end
end

#############################################################################
### TEST 5: MONTE CARLO CONVERGENCE RATE ####################################
#############################################################################

@testset "Convergence rate" begin
    Ndim = 2
    eps6 = ones(6)

    ref   = solvedVF(CUBE_POINTS, CUBE_FACES, Ndim, CUBE_QIN, CUBE_TIN, eps6)
    T_ref = collectT(ref)
    free  = findall(t -> t < 0, collectTin(ref))

    errs = Float64[]
    for rpe in (10_000, 1_000_000)
        d = solvedMC(CUBE_POINTS, CUBE_FACES, Ndim, CUBE_QIN, CUBE_TIN, eps6, rpe)
        @test abs(d.energy_error) < 1e-10
        push!(errs, sqrt(mean(abs2, collectT(d)[free] .- T_ref[free])))
    end

    ratio = errs[1] / errs[2]
    @test 6.5 < ratio < 15.0
end

#############################################################################
### TEST 6: OCCLUSION ON A NON-CONVEX ENCLOSURE #############################
#############################################################################

@testset "Occlusion in the L-duct" begin
    Ndim = 3
    d = tracedDomain(L_POINTS, L_FACES, Ndim, L_QIN, L_TIN, ones(14), 200_000)

    r = faceRanges(d)
    block = d.F_raw[r[9], r[12]]      # x = 2 wall against y = 2 wall

    # the segment between these walls leaves the solid at the reentrant corner,
    # so they are mutually hidden and the block must be structurally empty
    @test count(!iszero, block) == 0

    # sanity: faces 7 and 9 are both in the straight x-arm and mutually visible,
    # so a dead tracer returning an all-zero F would not pass the block above
    @test count(!iszero, d.F_raw[r[7], r[9]]) > 0
end