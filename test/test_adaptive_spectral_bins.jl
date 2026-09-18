# adaptiveSpectralBins tests
#
# 1. The three prototype spectra reach the requested tolerance with bin counts
#    in the prototype's range, and the model validates.
# 2. Pieces cover the sampled range without gaps; every piece's samples carry
#    κ inside the bin's interval.
# 3. Argument checking.
# 4. Isothermal enclosure with an adaptive model on the sigmoid spectrum.

println("\n" * "-"^60)
println("Testing adaptive spectral bins")
println("-"^60)

using Test
using RayTraceHeatTransfer
using RayTraceHeatTransfer: validate, n_pieces
using GeometryBasics, StaticArrays
using Random

κ_sigmoid(λ) = 0.01 + (100.0 - 0.01) / (1 + (4e-6 / λ)^6)

const AB_BANDS = [(2.7e-6, 0.05, 300.0), (4.3e-6, 0.03, 3000.0), (6.3e-6, 0.15, 40.0),
                  (15e-6, 0.20, 600.0), (1.9e-6, 0.03, 8.0)]
function κ_bands(λ)
    s = 0.05
    for (c, hw, p) in AB_BANDS
        x = log10(λ / c) / hw
        s += p / (1 + x^2)^2
    end
    return s
end

const AB_LINES = let rng = MersenneTwister(1)
    centres = 10 .^ (log10(1.5e-6) .+ (log10(30e-6) - log10(1.5e-6)) .* rand(rng, 400))
    peaks   = 10 .^ (log10(0.1) .+ 5.0 .* rand(rng, 400))
    widths  = 1e-4 .+ 2e-4 .* rand(rng, 400)
    collect(zip(centres, widths, peaks))
end
function κ_lines(λ)
    s = 1e-3
    for (c, hw, p) in AB_LINES
        x = log10(λ / c) / hw
        s += p / (1 + x^2)
    end
    return s
end

const AB_λ = 10 .^ range(log10(0.5e-6), log10(50e-6), length = 200001)
const AB_L = (0.01, 10.0)
const AB_T = (600.0, 2000.0)

@testset "tolerance reached, bin counts in range" begin
    expected = Dict("sigmoid"  => (κ_sigmoid.(AB_λ), (20, 60), (60, 200)),
                    "bands"    => (κ_bands.(AB_λ),   (20, 70), (60, 200)),
                    "lines"    => (κ_lines.(AB_λ),   (20, 80), (80, 300)))
    for (name, (κ, range3, range4)) in expected
        m3 = adaptiveSpectralBins(AB_λ, κ; tol = 1e-3, L_range = AB_L, T_range = AB_T)
        m4 = adaptiveSpectralBins(AB_λ, κ; tol = 1e-4, L_range = AB_L, T_range = AB_T)
        validate(m3, length(m3.κ_ref))
        validate(m4, length(m4.κ_ref))
        @test m3.achieved_error <= 1e-3
        @test m4.achieved_error <= 1e-4
        @test range3[1] <= length(m3.κ_ref) <= range3[2]
        @test range4[1] <= length(m4.κ_ref) <= range4[2]
        @test length(m4.κ_ref) > length(m3.κ_ref)
        @test issorted(m4.κ_lo)
        @test all(m4.κ_lo .<= m4.κ_ref .<= m4.κ_hi)
    end
end

@testset "pieces cover the range and respect bin intervals" begin
    κ = κ_bands.(AB_λ)
    m = adaptiveSpectralBins(AB_λ, κ; tol = 1e-3, L_range = AB_L, T_range = AB_T)
    @test m.edges[1] == AB_λ[1]
    @test m.edges[end] == AB_λ[end]
    @test all(diff(m.edges) .> 0)
    # every sample strictly inside a piece has κ in that piece's bin interval
    # (samples at a crossing edge may sit on either side by interpolation)
    bad = 0
    for p in 1:n_pieces(m)
        k = m.piece_bin[p]
        lo, hi = m.edges[p], m.edges[p+1]
        for i in searchsortedfirst(AB_λ, lo):searchsortedlast(AB_λ, hi)
            (AB_λ[i] == lo || AB_λ[i] == hi) && continue
            (m.κ_lo[k] <= κ[i] < m.κ_hi[k]) || (bad += 1)
        end
    end
    @test bad == 0
end

@testset "argument checking" begin
    λ = AB_λ[1:200:end]; κ = κ_sigmoid.(λ)
    @test_throws ArgumentError adaptiveSpectralBins(λ, κ[1:end-1]; tol = 1e-3, L_range = AB_L, T_range = AB_T)
    @test_throws ArgumentError adaptiveSpectralBins(reverse(λ), κ; tol = 1e-3, L_range = AB_L, T_range = AB_T)
    @test_throws ArgumentError adaptiveSpectralBins(λ, κ; tol = 0.0, L_range = AB_L, T_range = AB_T)
    @test_throws ArgumentError adaptiveSpectralBins(λ, κ; tol = 1e-3, L_range = (1.0, 0.1), T_range = AB_T)
    @test_throws ErrorException adaptiveSpectralBins(λ, κ; tol = 1e-3, L_range = AB_L, T_range = AB_T, max_bins = 3)
end

@testset "isothermal enclosure with an adaptive model" begin
    T_wall = 1000.0
    κ = κ_sigmoid.(AB_λ)
    model = adaptiveSpectralBins(AB_λ, κ; tol = 1e-2, L_range = (0.2, 1.5), T_range = (900.0, 1100.0))
    K = length(model.κ_ref)
    vertices = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    face = PolyVolume2D{Float64}(vertices, SVector(true, true, true, true), K, 1.0, 0.0)
    face.kappa_g   = copy(model.κ_ref)
    face.sigma_s_g = zeros(K)
    face.epsilon   = [fill(0.8, K) for _ in 1:4]
    face.T_in_w    = fill(T_wall, 4)
    face.q_in_w    = zeros(4)
    face.T_in_g    = -1.0
    face.q_in_g    = 0.0
    mesh = RayTracingDomain2D([face], [(5, 5)], verbose = false)
    mesh.spectral_model = model
    mesh(400_000; method = :exchange, verbose = false)
    smooth!(mesh; verbose = false)
    solveEquilibrium!(mesh, mesh.F_smooth; max_iters = 10_000, convergence_tol = 1e-14, verbose = false)
    err = maximum(abs(ff.T_g - T_wall) for ff in mesh.fine_mesh[1])
    @test err < 1e-9
end

println("✓ Adaptive spectral bins tests complete")