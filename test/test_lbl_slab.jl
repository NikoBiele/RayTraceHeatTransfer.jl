println("\n" * "-"^60)
println("Testing line-by-line spectral slab")
println("-"^60)

# Line-by-line plane-parallel slab reference, and the package's cavity-as-slab
# emulation of it — grey and with a synthetic line spectrum through adaptive
# spectral binning and one pathlength trace.
#
# Reference: absorbing–emitting, non-scattering gas between two black plates,
# temperature piecewise constant on Nx cells, exact exponential-integral
# quadrature per wavelength, Newton on the cell temperatures.

using Test
using RayTraceHeatTransfer
using RayTraceHeatTransfer: n_pieces
using GeometryBasics, StaticArrays
using SpecialFunctions      # expint
using LinearAlgebra
using Random
using Base.Threads

include("../examples/lbl_slab_reference.jl")

# ---- cavity-as-slab helpers -------------------------------------------------

# Wide cavity W×1 with black bottom (T1) and top (T2) plates and adiabatic,
# nearly reflecting sides; nx_h × nx_v volume cells. `kappa` is either a
# scalar (grey) or a per-bin vector (spectral).
function lbl_slab_domain(W, nx_h, nx_v, T1, T2, kappa)
    K = kappa isa AbstractVector ? length(kappa) : 1
    verts = SVector(Point2(0.0, 0.0), Point2(W, 0.0), Point2(W, 1.0), Point2(0.0, 1.0))
    if K == 1
        face = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), 1, kappa, 0.0)
        face.epsilon = [1.0, 0.01, 1.0, 0.01]   # sides: adiabatic, nearly reflecting (ε = 0 is undefined)
    else
        face = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), K, 1.0, 0.0)
        face.kappa_g   = copy(kappa)
        face.sigma_s_g = zeros(K)
        face.epsilon   = [fill(1.0, K), fill(0.01, K), fill(1.0, K), fill(0.01, K)]
    end
    face.T_in_w = [T1, -1.0, T2, -1.0]
    face.q_in_w = zeros(4)
    face.T_in_g = -1.0
    face.q_in_g = 0.0
    return RayTracingDomain2D([face], [(nx_h, nx_v)], verbose = false)
end

# Volume cells of the centre column, bottom to top.
function lbl_centre_column(mesh, W, nx_h)
    ic  = (nx_h + 1) ÷ 2
    col = [f for f in mesh.fine_mesh[1] if abs(f.midPoint[1] - (ic - 0.5) * W / nx_h) < 1e-9]
    sort!(col, by = f -> f.midPoint[2])
    return col
end

# (ψ, T-profile) of the centre column; ψ from the bottom-wall element's net source.
function lbl_column_result(mesh, W, nx_h, T1, T2)
    col   = lbl_centre_column(mesh, W, nx_h)
    T_pkg = [f.T_g for f in col]
    f0    = col[1]
    ψ_pkg = (f0.q_w[1] / f0.area[1]) / (LBL_σ * (T1^4 - T2^4))
    return ψ_pkg, T_pkg
end

# ---- tests --------------------------------------------------------------------

@testset "LBL slab reference vs Heaslet & Warming (grey)" begin
    # ψ_b(τ_L) from Heaslet & Warming (1965), Table 13.1 in Modest & Mazumder,
    # Radiative Heat Transfer, 4th ed. (n = 1)
    HW = [(0.1, 0.9157), (0.2, 0.8491), (0.3, 0.7934), (0.4, 0.7458),
          (0.5, 0.7040), (0.6, 0.6672), (0.8, 0.6046), (1.0, 0.5532),
          (1.5, 0.4572), (2.0, 0.3900), (2.5, 0.3401), (3.0, 0.3016),
          (5.0, 0.2077)]

    T1, T2 = 1000.0, 500.0
    λg = 10 .^ range(log10(1e-7), log10(1e-3), length = 801)
    for (τL, ψ_ref) in HW
        _, q, it = lbl_slab_equilibrium(λg, fill(τL, length(λg)), 1.0, T1, T2; Nx = 100, tol = 1e-9)
        ψ = q / (LBL_σ * (T1^4 - T2^4))
        @test abs(ψ - ψ_ref) < 1e-3
        @test it < 30
    end
    # optically thick limit: ψ_b ≈ (4/3) / (1.42089 + τ_L) for τ_L ≫ 1 (same table)
    τL = 10.0
    _, q, _ = lbl_slab_equilibrium(λg, fill(τL, length(λg)), 1.0, T1, T2; Nx = 200, tol = 1e-9)
    @test abs(q / (LBL_σ * (T1^4 - T2^4)) - (4 / 3) / (1.42089 + τL)) < 2e-3
end

@testset "Cavity-as-slab emulation, grey gas" begin
    T1, T2 = 1000.0, 500.0
    W, NX  = 20.0, 16
    nx_h   = round(Int, W * NX / 5)              # 5:1 cell aspect
    τL     = 1.0
    mesh   = lbl_slab_domain(W, nx_h, NX, T1, T2, τL)
    mesh(10_000_000; method = :pathlength, verbose = false)
    smooth!(mesh; verbose = false)
    solveEquilibrium!(mesh, mesh.F_smooth; max_iters = 20_000, convergence_tol = 1e-12, verbose = false)
    ψ_pkg, T_pkg = lbl_column_result(mesh, W, nx_h, T1, T2)

    λg = 10 .^ range(log10(1e-7), log10(1e-3), length = 801)
    T_lbl, q_lbl, _ = lbl_slab_equilibrium(λg, fill(τL, length(λg)), 1.0, T1, T2; Nx = NX, tol = 1e-9)
    ψ_lbl = q_lbl / (LBL_σ * (T1^4 - T2^4))

    @test abs(ψ_pkg - ψ_lbl) < 1e-2
    @test maximum(abs.(T_pkg .- T_lbl)) < 5.0
end

@testset "Cavity-as-slab emulation, synthetic line spectrum + adaptive bins" begin
    ψ_TOL   = 3e-3      # |ψ_pkg − ψ_LBL|
    ΔT_TOL  = 8.0       # max |T_pkg − T_LBL| [K]
    MAX_BIN = 40        # adaptive bins at tol = 1e-2 (20 on the 200k-point grid)
    N_RAYS  = 2_000_000

    T1, T2 = 1000.0, 500.0
    NX, W, NX_H = 32, 1000.0, 5

    # 400 Lorentzian lines on a weak continuum, strengths over five decades
    λ = 10 .^ range(log10(1e-8), log10(1e-2), length = 60_001)
    lines = let rng = MersenneTwister(1)
        centres = 10 .^ (log10(1.5e-6) .+ (log10(30e-6) - log10(1.5e-6)) .* rand(rng, 400))
        peaks   = 10 .^ (log10(0.1) .+ 5.0 .* rand(rng, 400))
        widths  = 1e-4 .+ 2e-4 .* rand(rng, 400)
        collect(zip(centres, widths, peaks))
    end
    κ = [1e-3 + sum(p / (1 + (log10(x / c) / hw)^2) for (c, hw, p) in lines) for x in λ]
    κ .*= 0.1

    T_lbl, q_lbl, it = lbl_slab_equilibrium(λ, κ, 1.0, T1, T2; Nx = NX, tol = 1e-8)
    ψ_lbl = q_lbl / (LBL_σ * (T1^4 - T2^4))
    @test it < 30

    model = adaptiveSpectralBins(λ, κ; tol = 1e-2, L_range = (1 / NX, 3.0), T_range = (T2, T1))
    K = length(model.κ_ref)
    @test K <= MAX_BIN
    @test model.achieved_error <= 1e-2
    @test n_pieces(model) > K                    # bins are unions of λ-pieces

    mesh = lbl_slab_domain(W, NX_H, NX, T1, T2, model.κ_ref)
    mesh(N_RAYS; method = :pathlength, verbose = false)
    mesh.spectral_model = model
    smooth!(mesh; verbose = false, k_dykstra = 400)
    solveEquilibrium!(mesh, mesh.F_smooth; max_iters = 20_000, convergence_tol = 1e-12, verbose = false)
    ψ_pkg, T_pkg = lbl_column_result(mesh, W, NX_H, T1, T2)

    @test abs(ψ_pkg - ψ_lbl) < ψ_TOL
    @test maximum(abs.(T_pkg .- T_lbl)) < ΔT_TOL
    @test all(T2 .< T_pkg .< T1)                 # profile bracketed by the plates
end

println("✓ Line-by-line spectral slab tests complete")