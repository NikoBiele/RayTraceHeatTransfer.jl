# Arbitrary precision on the 3D analytical view factor path: the exact Clausen
# function, the parallel-edge term, σ in any type, and a cube solved end to end
# in BigFloat, where reciprocity, energy conservation and temperatures all
# follow the working precision.

println("\n" * "-"^60)
println("Testing arbitrary precision (3D view factors and grey solve)")
println("-"^60)

using Test
using RayTraceHeatTransfer
using Random
const RTHT = RayTraceHeatTransfer

# Unit cube: hot face 1000 K, cold face 0 K, four sides in radiative equilibrium
cube_domain(T, Ndim) = ViewFactorDomain3D(
    T.([0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1]),   # vertices in type T
    [1 2 4 3; 5 6 8 7; 1 5 7 3; 2 6 8 4; 3 4 8 7; 1 2 6 5],          # six quad faces
    Ndim, zeros(T, 6), T.([1000, 0, -1, -1, -1, -1]), ones(T, 6))     # q, T (−1 = unknown), ε

# View factors and grey solve on F_raw, all in the domain's element type
function solved_cube(T, Ndim)
    d = cube_domain(T, Ndim)
    d(; verbose = false)
    solveEquilibrium!(d, d.F_raw; verbose = false)
    return d
end

# Relative reciprocity defect max|AᵢFᵢⱼ − AⱼFⱼᵢ| / max(AᵢFᵢⱼ) of F_raw
function reciprocity_defect(d)
    A  = [sf.area for f in d.facesMesh for sf in f.subFaces]
    AF = A .* d.F_raw
    return maximum(abs.(AF .- AF')) / maximum(AF)
end

@testset "Clausen function" begin
    setprecision(BigFloat, 128) do
        ϵ = eps(BigFloat)
        # Cl₂(π/2) is Catalan's constant, Cl₂(π) = 0
        @test abs(RTHT.Cl3D(big(pi) / 2) - BigFloat(Base.MathConstants.catalan)) < 10ϵ
        @test abs(RTHT.Cl3D(big(pi))) < 10ϵ
        # duplication formula Cl₂(2θ) = 2Cl₂(θ) − 2Cl₂(π − θ)
        θ = big"0.7"
        @test abs(RTHT.Cl3D(2θ) - 2RTHT.Cl3D(θ) + 2RTHT.Cl3D(big(pi) - θ)) < 10ϵ
        # Float64 agrees with the 128-bit value over a full period
        θs = range(1e-6, 2π - 1e-6; length = 2001)
        @test maximum(abs(RTHT.Cl3D(t) - Float64(RTHT.Cl3D(big(t)))) for t in θs) < 1e-14
    end
end

@testset "Parallel-edge term" begin
    setprecision(BigFloat, 128) do
        rng  = MersenneTwister(1)
        args = [(randn(rng), randn(rng), abs(randn(rng))) for _ in 1:1000]           # general (s, l, d)
        near = [(2 + rand(rng), rand(rng), 1e-3 * rand(rng)) for _ in 1:1000]      # old clamp region
        for set in (args, near)
            err = maximum(abs(RTHT.fParallel3D(s, l, d, 0.0) -
                              Float64(RTHT.fParallel3D(big(s), big(l), big(d), big(0)))) for (s, l, d) in set)
            @test err < 1e-13
        end
    end
end

@testset "Stefan–Boltzmann constant" begin
    setprecision(BigFloat, 256) do
        @test Float64(RTHT.stefan_boltzmann(BigFloat)) == RTHT.STEFAN_BOLTZMANN   # correctly rounded literal
        @test RTHT.stefan_boltzmann(BigFloat) isa BigFloat
    end
end

@testset "Cube end to end in BigFloat" begin
    d64 = solved_cube(Float64, 2)
    setprecision(BigFloat, 128) do
        ϵ = eps(BigFloat)
        d = solved_cube(BigFloat, 2)
        # precision is carried through, not truncated anywhere
        @test d.F_raw isa Matrix{BigFloat}
        @test d.facesMesh[1].subFaces[1].vertices[1] isa RTHT.Point3{BigFloat}
        @test d.facesMesh[3].subFaces[1].T_w isa BigFloat
        # analytical reciprocity survives row normalisation to working precision
        @test reciprocity_defect(d) < 100ϵ
        # energy conservation to working precision
        @test abs(d.energy_error) < 100ϵ
        # Float64 temperatures agree with the 128-bit ones to Float64 round-off
        T64  = [sf.T_w for f in d64.facesMesh for sf in f.subFaces]
        Tbig = [sf.T_w for f in d.facesMesh   for sf in f.subFaces]
        @test maximum(abs.(T64 .- Tbig)) < 1e-9
        # the summary prints, with a short energy error
        s = sprint(show, MIME"text/plain"(), d)
        @test occursin("relative conservation error", s)
        @test length(s) < 2000
    end
end

println("✓ Precision tests complete")