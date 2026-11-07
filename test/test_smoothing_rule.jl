println("\n" * "-"^60)
println("Testing the smoothing cost rule (Dykstra vs alternating projections)")
println("-"^60)

using Test
using RayTraceHeatTransfer
using RayTraceHeatTransfer: predict_rho_data, get_w, DualSolver, dykstra_rate_bound, DkAP
using GeometryBasics, StaticArrays
using LinearAlgebra

# Rectangular black cavity Wc x H with a grey gas; the two long side walls are
# the blind set, so Wc/H -> 0 with kappa*Wc -> 0 approaches the parallel-plate limit.
function rule_cavity(Wc, H, nx_h, nx_v, kappa)
    verts = SVector(Point2(0.0, 0.0), Point2(Wc, 0.0), Point2(Wc, H), Point2(0.0, H))
    face  = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), 1, kappa, 0.0)
    face.epsilon = [1.0, 1.0, 1.0, 1.0]
    face.T_in_w  = [1000.0, 500.0, 500.0, 500.0]
    face.q_in_w  = zeros(4)
    face.T_in_g  = -1.0
    face.q_in_g  = 0.0
    return RayTracingDomain2D([face], [(nx_h, nx_v)], verbose = false)
end

@testset "Dykstra rate bound" begin
    println("  rate bound on a dense and a thin-slit support")
    w = [1.0, 2.0, 3.0, 4.0]
    dual = DualSolver(w)
    @test dykstra_rate_bound(ones(4, 4), dual) == 0.0          # full support: no bound
    @test dykstra_rate_bound(Matrix(1.0I, 4, 4), dual) > 0.5   # diagonal only: slow

    mesh = rule_cavity(0.002, 10.0, 2, 40, 0.05)
    mesh(2_000_000; verbose = false)
    w  = get_w(mesh); w ./= minimum(w)
    rho_d = dykstra_rate_bound(Matrix(mesh.F_raw), DualSolver(w))
    @test 0.9 < rho_d < 1.0                                    # banded support, ~1 - 1/nx_v
end

@testset "Branch 3: thin slit runs Dykstra to its floor" begin
    println("  thin slit (near parallel plates), sparse input densified, Dykstra expected")
    mesh = rule_cavity(0.002, 10.0, 2, 40, 0.05)
    mesh(2_000_000; verbose = false)
    pred = predict_rho_data(mesh)
    @test pred.rho_est > 0.99
    stats = smooth!(mesh; verbose = false)
    @test stats.k_dykstra[1] > 1                               # k_dykstra = 1 was overridden
    @test stats.k_pcg_tot[1] > stats.k_pcg_max[1]              # more than one round ran
    @test stats.converged[1]
    @test stats.delta_smooth[1] < 1e-12
    @test all(isfinite, mesh.F_smooth)
    @test minimum(mesh.F_smooth) >= 0
    @test maximum(abs.(sum(mesh.F_smooth, dims = 2) .- 1)) < 1e-12
end

@testset "Opt-out: k_dykstra = 0 gives pure alternating projections" begin
    println("  thin slit with k_dykstra = 0 and a small k_ap cap")
    mesh = rule_cavity(0.002, 10.0, 2, 40, 0.05)
    mesh(2_000_000; verbose = false)
    stats = @test_logs (:warn,) match_mode=:any smooth!(mesh; k_dykstra = 0, k_ap = 50, verbose = false)
    @test stats.k_dykstra[1] == 0
    @test stats.k_pcg_tot[1] == 0
    @test stats.k_ap[1] == 50
    @test !stats.converged[1]                                  # too slow for 50 steps, as predicted
end

@testset "Ordinary cavity: the branch taken agrees with the rule" begin
    println("  ordinary cavity, dense, branch checked against exchange rate and cost ratio")
    mesh = rule_cavity(2.0, 3.0, 4, 6, 0.05)
    mesh(2_000_000; verbose = false)
    pred = predict_rho_data(mesh)
    @test pred.rho_est < 0.95                                  # far from the bipartite limit
    w  = get_w(mesh); w ./= minimum(w)
    N  = length(w)
    rho_d = dykstra_rate_bound(Matrix(mesh.F_raw), DualSolver(w))
    nnz_raw = count(!iszero, mesh.F_raw)
    exchange = log(rho_d) / log(pred.rho_est)
    stats = smooth!(mesh; verbose = false)
    @test stats.converged[1]
    @test stats.delta_smooth[1] < 1e-12
    @test stats.k_dykstra[1] >= 1
    # PCG counts are flat across rounds, so k_pcg_max bounds the first-round
    # count from above and the cost ratio built from it bounds the true one.
    cost_ratio = (5 + stats.k_pcg_max[1]) * N^2 / (3 * nnz_raw)
    if exchange > cost_ratio
        @test stats.k_dykstra[1] > 1                           # Dykstra had to pay
    end
    if stats.k_dykstra[1] == 1
        @test exchange <= cost_ratio                           # a single round means it did not
    end
end

@testset "Both branches forced on one matrix through rho_ap" begin
    println("  DkAP called directly with a slow and a fast a priori AP rate")
    mesh = rule_cavity(2.0, 3.0, 4, 6, 0.05)
    mesh(2_000_000; verbose = false)
    w  = get_w(mesh); w ./= minimum(w)
    F  = Matrix(mesh.F_raw)
    n_s = length(mesh.surface_mapping)
    nz_over_N = count(!iszero, F) / length(w)
    # rho_ap far below the Dykstra bound: a round never pays, prescribed rounds kept
    out2 = DkAP(F, w, n_s; k_dykstra = 3, k_ap = 1_000, verbose = false, nz_over_N = nz_over_N, rho_ap = 0.3)
    @test out2[5] == 3                                         # k_dykstra honoured
    @test out2[2]                                              # converged
    # rho_ap near 1: every round pays, Dykstra runs to its floor
    out3 = DkAP(F, w, n_s; k_dykstra = 3, k_ap = 1_000, verbose = false, nz_over_N = nz_over_N, rho_ap = 0.9999)
    @test out3[5] > 3                                          # prescribed count overridden
    @test out3[2]
    @test out3[7] > out2[7]                                    # more PCG work in total
end

println("✓ Smoothing cost rule tests passed")