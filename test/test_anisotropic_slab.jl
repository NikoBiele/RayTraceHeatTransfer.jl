# Anisotropic slab: deterministic 1D discrete-ordinates reference, checked against
# Heaslet & Warming and the grey line-by-line reference, then used to validate the
# grey directional solver (Henyey–Greenstein scattering) on a cavity emulating the slab.

println("\n" * "-"^60)
println("Testing the directional solver against the anisotropic slab reference")
println("-"^60)

using Test
using RayTraceHeatTransfer
using GeometryBasics, StaticArrays, LinearAlgebra

include("../examples/anisotropic_slab_reference.jl")

# Cavity emulating the slab: bottom wall T1, top wall T2, black; side walls adiabatic.
# With mirror_sides the side walls are specular, which makes a cavity of any width equivalent
# to the infinite slab by symmetry (ε_side > 0 is needed for a radiative-equilibrium surface,
# so that fraction of their interaction stays diffuse); otherwise they reflect diffusely and
# the cavity must be very wide.
function aslab_domain(W, nx_h, nx_v, T1, T2, κ, σ_s, g; ε_side = 0.01, mirror_sides = true,
                      bins = AngularBins(16, 4))
    verts = SVector(Point2(0.0, 0.0), Point2(W, 0.0), Point2(W, 1.0), Point2(0.0, 1.0))
    face  = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), 1, Float64(κ), Float64(σ_s))
    face.epsilon = [1.0, ε_side, 1.0, ε_side]
    face.T_in_w  = [T1, -1.0, T2, -1.0]
    face.q_in_w  = zeros(4)
    face.T_in_g  = -1.0
    face.q_in_g  = 0.0
    face.phase   = g == 0 ? IsotropicScattering() : HenyeyGreenstein(g)
    side = mirror_sides ? SpecularReflection(1.0) : DiffuseReflection()
    face.reflection = [DiffuseReflection(), side, DiffuseReflection(), side]
    mesh = RayTracingDomain2D([face], [(nx_h, nx_v)], verbose = false)
    mesh.directional_model = bins
    return mesh
end

# row-averaged gas temperatures (bottom to top) and ψ from the bottom-wall elements
function aslab_result(mesh, nx_v, T1, T2)
    T = zeros(nx_v); cnt = zeros(Int, nx_v)
    for f in mesh.fine_mesh[1]
        r = clamp(floor(Int, f.midPoint[2] * nx_v) + 1, 1, nx_v)
        T[r] += f.T_g; cnt[r] += 1
    end
    T ./= cnt
    q = 0.0; len = 0.0
    for ((ci, fi, wi), _) in mesh.surface_mapping
        f  = mesh.fine_mesh[ci][fi]
        p1 = f.vertices[wi]; p2 = f.vertices[mod1(wi + 1, length(f.vertices))]
        if p1[2] < 1e-9 && p2[2] < 1e-9                       # bottom wall
            q += f.q_w[wi]; len += f.area[wi]
        end
    end
    return T, (q / len) / (ASLAB_σ * (T1^4 - T2^4))
end

@testset "reference vs Heaslet & Warming (g = 0)" begin
    T1, T2 = 1000.0, 500.0
    for (τL, ψ_ref) in [(0.5, 0.7040), (1.0, 0.5532), (2.0, 0.3900), (10.0, 0.1167)]
        _, q = anisotropic_slab_equilibrium(τL, 0.0, 0.0, 1.0, T1, T2; Nx = 200, n_half = 8)
        @test abs(q / (ASLAB_σ * (T1^4 - T2^4)) - ψ_ref) < 5e-4
    end
end

@testset "reference: invariances" begin
    T1, T2 = 1000.0, 500.0
    # in radiative equilibrium isotropic scattering is indistinguishable from absorption
    Ta, qa = anisotropic_slab_equilibrium(1.0, 0.0, 0.0, 1.0, T1, T2; Nx = 32, n_half = 16)
    Ts, qs = anisotropic_slab_equilibrium(0.2, 0.8, 0.0, 1.0, T1, T2; Nx = 32, n_half = 16)
    @test maximum(abs.(Ta .- Ts)) < 1e-8 && abs(qa - qs) / qa < 1e-12
    # isothermal walls, grey mirror walls, forward scattering: nothing happens
    Ti, qi = anisotropic_slab_equilibrium(0.2, 0.8, 0.8, 1.0, 800.0, 800.0; ε1 = 0.5, ε2 = 0.5,
                                          specularity = 1.0, Nx = 16, n_half = 12)
    @test maximum(abs.(Ti .- 800.0)) < 1e-7 && abs(qi) < 1e-6
    # forward scattering raises the flux and flattens the profile
    Tg, qg = anisotropic_slab_equilibrium(0.2, 0.8, 0.8, 1.0, T1, T2; Nx = 32, n_half = 32)
    @test qg > 1.2 * qs
    @test Tg[1] < Ts[1] - 10 && Tg[end] > Ts[end] + 10
    # converged in the number of ordinates
    Th, qh = anisotropic_slab_equilibrium(0.2, 0.8, 0.8, 1.0, T1, T2; Nx = 32, n_half = 16)
    @test maximum(abs.(Th .- Tg)) < 1e-2 && abs(qh - qg) / qg < 1e-5
end

@testset "reference vs the grey line-by-line reference (same cells)" begin
    T1, T2, NX = 1000.0, 500.0, 32
    λg = 10 .^ range(log10(1e-7), log10(1e-3), length = 801)
    T_lbl, q_lbl, _ = lbl_slab_equilibrium(λg, fill(1.0, length(λg)), 1.0, T1, T2; Nx = NX, tol = 1e-9)
    T_do, q_do = anisotropic_slab_equilibrium(1.0, 0.0, 0.0, 1.0, T1, T2; Nx = NX, n_half = 32)
    @test maximum(abs.(T_do .- T_lbl)) < 0.1
    @test abs(q_do - q_lbl) / q_lbl < 5e-4
end

@testset "grey directional solver vs the reference (Henyey–Greenstein, g = 0.8)" begin
    T1, T2 = 1000.0, 500.0
    κ, σ_s, g = 0.2, 0.8, 0.8
    W, NX_H, NX = 2.0, 4, 32
    T_ref, q_ref = anisotropic_slab_equilibrium(κ, σ_s, g, 1.0, T1, T2; Nx = NX, n_half = 32)
    T_iso, q_iso = anisotropic_slab_equilibrium(κ, σ_s, 0.0, 1.0, T1, T2; Nx = NX, n_half = 32)
    ψ_ref = q_ref / (ASLAB_σ * (T1^4 - T2^4))
    ψ_iso = q_iso / (ASLAB_σ * (T1^4 - T2^4))

    mesh = aslab_domain(W, NX_H, NX, T1, T2, κ, σ_s, g)
    mesh(4_000_000; method = :pathlength, verbose = false)
    stats = smooth!(mesh; verbose = false)
    @test all(stats.converged)
    solveEquilibrium!(mesh, mesh.F_smooth; verbose = false)
    T_pkg, ψ_pkg = aslab_result(mesh, NX, T1, T2)

    ΔT = maximum(abs.(T_pkg .- T_ref))
    @test abs(mesh.energy_error) < 1e-10
    @test ΔT < 3.0
    @test abs(ψ_pkg - ψ_ref) / ψ_ref < 0.02
    # the solution is the anisotropic one, not the isotropic one
    @test maximum(abs.(T_pkg .- T_ref)) < maximum(abs.(T_pkg .- T_iso)) / 3
    @test abs(ψ_pkg - ψ_ref) < abs(ψ_pkg - ψ_iso) / 3
end

println("✓ Directional solver against reference solution tests complete")