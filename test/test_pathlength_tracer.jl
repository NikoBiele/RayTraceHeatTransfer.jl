# Pathlength tracer: one geometric trace, exchange factors for any κ.

println("\n" * "-"^60)
println("Testing method = :pathlength and exchangeFactors!")
println("-"^60)

using Test
using RayTraceHeatTransfer
using GeometryBasics, StaticArrays
using LinearAlgebra, SparseArrays

function _pl_face(κ, ε_wall, T_wall)
    n = length(κ)
    verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    face = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), n, 1.0, 0.0)
    face.kappa_g   = copy(κ)
    face.sigma_s_g = zeros(n)
    face.epsilon   = [fill(ε_wall, n) for _ in 1:4]
    face.T_in_w    = fill(T_wall, 4)
    face.q_in_w    = zeros(4)
    face.T_in_g    = -1.0
    face.q_in_g    = 0.0
    return face
end

const PL_κ = [0.05, 3.0, 0.3, 10.0]

@testset "store geometry" begin
    mesh = RayTracingDomain2D([_pl_face(PL_κ, 1.0, 1000.0)], [(4, 4)], verbose = false)
    mesh(200_000; method = :pathlength, verbose = false)
    store = mesh.path_store
    ns, nv = store.num_surfaces, store.num_volumes
    @test all(ns .< store.seg_cell .<= ns + nv)
    @test all(1 .<= store.ray_end .<= ns)
    nr = length(store.ray_emitter)
    ray_len = [sum(store.seg_len[store.ray_start[r]:store.ray_start[r+1]-1]) for r in 1:nr]
    inplane = ray_len .* [hypot(d[1], d[2]) for d in store.ray_dir]
    @test maximum(inplane) <= sqrt(2) * (1 + 1e-5)
    vol = store.ray_emitter .> ns
    @test all(store.seg_cell[store.ray_start[1:end-1]][vol] .== store.ray_emitter[vol])
end

@testset "agreement with the exchange tracer" begin
    mesh_p = RayTracingDomain2D([_pl_face(PL_κ, 1.0, 1000.0)], [(4, 4)], verbose = false)
    mesh_p(2_000_000; method = :pathlength, verbose = false, seeds = 100)
    exchangeFactors!(mesh_p; verbose = false)
    mesh_e = RayTracingDomain2D([_pl_face(PL_κ, 1.0, 1000.0)], [(4, 4)], verbose = false)
    mesh_e(2_000_000; method = :exchange, verbose = false, seeds = 100)
    for k in 1:4
        Fp, Fe = mesh_p.F_raw[k], mesh_e.F_raw[k]
        @test maximum(abs.(sum(Fp, dims = 2) .- 1)) < 1e-12
        @test norm(Fp - Fe) / norm(Fe) < 2e-2
    end
end

@testset "isothermal enclosure, exact" begin
    for (ε_wall, model) in ((1.0, PlanckBands([5e-7, 2e-6, 4e-6, 8e-6, 5e-5])),
                            (0.5, ConstantWeights([0.4, 0.3, 0.2, 0.1])))
        mesh = RayTracingDomain2D([_pl_face(PL_κ, ε_wall, 1000.0)], [(5, 5)], verbose = false)
        mesh.spectral_model = model
        mesh(1_000_000; method = :pathlength, verbose = false)
        exchangeFactors!(mesh; verbose = false)
        smooth!(mesh; verbose = false)
        solveEquilibrium!(mesh, mesh.F_smooth; max_iters = 10_000, convergence_tol = 1e-14, verbose = false)
        @test maximum(abs(ff.T_g - 1000.0) for ff in mesh.fine_mesh[1]) < 1e-9
    end
end

@testset "transparent bin and re-evaluation without retracing" begin
    mesh = RayTracingDomain2D([_pl_face(PL_κ, 0.8, 1000.0)], [(5, 5)], verbose = false)
    mesh.spectral_model = PlanckBands([5e-7, 2e-6, 4e-6, 8e-6, 5e-5])
    mesh(1_000_000; method = :pathlength, verbose = false)
    # change κ after tracing: bin 1 transparent, bin 4 stronger
    κ2 = [0.0, 3.0, 0.3, 30.0]
    for ff in mesh.fine_mesh[1]
        ff.kappa_g = copy(κ2)
    end
    exchangeFactors!(mesh; verbose = false)
    ns = length(mesh.surface_mapping)
    @test all(mesh.F_raw[1][:, ns+1:end] .== 0)
    smooth!(mesh; verbose = false)
    solveEquilibrium!(mesh, mesh.F_smooth; max_iters = 10_000, convergence_tol = 1e-14, verbose = false)
    @test maximum(abs(ff.T_g - 1000.0) for ff in mesh.fine_mesh[1]) < 1e-9
end

println("✓ Pathlength tracer and exchange factors tests complete")