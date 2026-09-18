# PiecewiseBands tests
#
# 1. One piece per band reproduces PlanckBands exactly (fractions and
#    derivatives), so the piecewise model is a strict generalisation.
# 2. Many pieces scattered over few bins: fractions sum to 1 and derivatives
#    sum to 0 at any temperature; derivatives match finite differences.
# 3. Validation rejects malformed models.
# 4. Isothermal enclosure with a PiecewiseBands model solves to T_wall.

println("\n" * "-"^60)
println("Testing piecewise bands spectral model")
println("-"^60)

using Test
using RayTraceHeatTransfer
using RayTraceHeatTransfer: fill_bin_fractions!, bin_fraction_derivatives!, validate
using GeometryBasics, StaticArrays
using Random

@testset "PiecewiseBands with one piece per band equals PlanckBands" begin
    limits = 10 .^ range(log10(1e-7), log10(1e-3), length = 21)
    K = length(limits) - 1
    m_pb = PlanckBands(limits)
    m_pw = PiecewiseBands(copy(limits), collect(1:K), ones(K), zeros(K), ones(K), 0.0)
    validate(m_pw, K)
    row_pb = zeros(K); row_pw = zeros(K)
    for T in (300.0, 1000.0, 3000.0)
        fill_bin_fractions!(row_pb, m_pb, T)
        fill_bin_fractions!(row_pw, m_pw, T)
        @test maximum(abs.(row_pb .- row_pw)) < 1e-11
        @test sum(row_pw) ≈ 1.0 atol = 1e-12
        bin_fraction_derivatives!(row_pb, m_pb, T)
        bin_fraction_derivatives!(row_pw, m_pw, T)
        @test maximum(abs.(row_pb .- row_pw)) < 1e-14 * maximum(abs.(row_pb))
    end
end

@testset "Scattered pieces: partition of unity and derivatives" begin
    rng = MersenneTwister(3)
    P, K = 400, 7
    edges = sort(10 .^ (log10(5e-7) .+ (log10(5e-5) - log10(5e-7)) .* rand(rng, P + 1)))
    piece_bin = rand(rng, 1:K, P)
    piece_bin[1:K] .= 1:K                       # every bin owns at least one piece
    m = PiecewiseBands(edges, piece_bin, ones(K), zeros(K), ones(K), 0.0)
    validate(m, K)
    row = zeros(K); drow = zeros(K); rp = zeros(K); rm = zeros(K)
    for T in (400.0, 1500.0, 5800.0)
        fill_bin_fractions!(row, m, T)
        @test sum(row) ≈ 1.0 atol = 1e-12
        @test all(row .>= 0)
        bin_fraction_derivatives!(drow, m, T)
        @test abs(sum(drow)) < 1e-12 * maximum(abs.(drow))
        h = 1e-3 * T
        fill_bin_fractions!(rp, m, T + h)
        fill_bin_fractions!(rm, m, T - h)
        fd = (rp .- rm) ./ (2h)
        @test maximum(abs.(fd .- drow)) < 1e-6 * maximum(abs.(drow))
    end
end

@testset "PiecewiseBands validation" begin
    edges = [1e-6, 2e-6, 3e-6, 4e-6]
    @test_throws ErrorException validate(PiecewiseBands(edges, [1, 2, 3], ones(3), zeros(3), ones(3), 0.0), 4)
    @test_throws ErrorException validate(PiecewiseBands([1e-6, 3e-6, 2e-6, 4e-6], [1, 2, 3], ones(3), zeros(3), ones(3), 0.0), 3)
    @test_throws ErrorException validate(PiecewiseBands(edges, [1, 1, 2], ones(3), zeros(3), ones(3), 0.0), 3)
    @test_throws ErrorException validate(PiecewiseBands(edges, [1, 2], ones(3), zeros(3), ones(3), 0.0), 3)
end

@testset "Isothermal enclosure with PiecewiseBands" begin
    T_wall = 1000.0
    κ_ref  = [0.5, 1.0, 1.5]
    # six pieces spanning the Planck-relevant range, folded onto three bins
    edges     = [5e-7, 1.5e-6, 2.5e-6, 4e-6, 7e-6, 1.5e-5, 5e-5]
    piece_bin = [1, 2, 3, 3, 2, 1]
    model = PiecewiseBands(edges, piece_bin, κ_ref, κ_ref ./ 1.2, κ_ref .* 1.2, 0.0)

    vertices = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    for ε_wall in (1.0, 0.5)
        face = PolyVolume2D{Float64}(vertices, SVector(true, true, true, true), 3, 1.0, 0.0)
        face.kappa_g   = copy(κ_ref)
        face.sigma_s_g = zeros(3)
        face.epsilon   = [fill(ε_wall, 3) for _ in 1:4]
        face.T_in_w    = fill(T_wall, 4)
        face.q_in_w    = zeros(4)
        face.T_in_g    = -1.0
        face.q_in_g    = 0.0
        mesh = RayTracingDomain2D([face], [(5, 5)], verbose = false)
        mesh.spectral_model = model
        mesh(1_000_000; method = :exchange, verbose = false)
        smooth!(mesh; verbose = false)
        solveEquilibrium!(mesh, mesh.F_smooth;
                          max_iters = 10_000, convergence_tol = 1e-14, verbose = false)
        err = maximum(abs(ff.T_g - T_wall) for ff in mesh.fine_mesh[1])
        @test err < 1e-9
    end
end

println("✓ Piecewise bands tests complete")