println("\n" * "-"^60)
println("Testing directional solvers")
println("-"^60)

using Test
using RayTraceHeatTransfer
using Random, LinearAlgebra, SparseArrays
using GeometryBasics, StaticArrays

@testset "grey directional solver" begin
    M = RayTraceHeatTransfer
    verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    dm = AngularBins(8, 2)
    function cavity(T_walls; phase = IsotropicScattering(), reflection = DiffuseReflection())
        f = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), 1, 0.2, 0.8)
        f.epsilon = fill(0.5, 4); f.T_in_w = copy(T_walls); f.q_in_w = zeros(4)
        f.T_in_g = -1.0; f.q_in_g = 0.0
        f.phase = phase
        f.reflection = fill(reflection, 4)
        mesh = RayTracingDomain2D([f], [(4, 4)], verbose = false)
        mesh.directional_model = dm
        mesh(400_000; method = :pathlength, seeds = 100, verbose = false)
        return mesh
    end
    gasT(mesh) = [mesh.fine_mesh[ci][fi].T_g for ((ci, fi), _) in sort(collect(mesh.volume_mapping), by = last)]
    gasj(mesh) = [mesh.fine_mesh[ci][fi].j_g for ((ci, fi), _) in sort(collect(mesh.volume_mapping), by = last)]

    # 1. default descriptors reproduce the grey solver exactly, on raw and on smoothed exchange factors
    mesh = cavity([1000.0, 500.0, 500.0, 500.0])
    M.equilibriumGrey2D!(mesh, Matrix(mesh.F_raw); verbose = false)
    T_ref, j_ref = gasT(mesh), gasj(mesh)
    solveEquilibrium!(mesh, mesh.F_raw; verbose = false)
    @test maximum(abs.(gasT(mesh) .- T_ref)) < 1e-7
    @test maximum(abs.(gasj(mesh) .- j_ref) ./ j_ref) < 1e-9
    @test abs(mesh.energy_error) < 1e-10

    @test_throws ErrorException solveEquilibrium!(mesh, mesh.F_smooth; verbose = false)   # not smoothed yet
    smooth!(mesh; verbose = false)
    M.equilibriumGrey2D!(mesh, mesh.F_smooth; verbose = false)
    T_diffuse = gasT(mesh)
    solveEquilibrium!(mesh, mesh.F_smooth; verbose = false)
    @test maximum(abs.(gasT(mesh) .- T_diffuse)) < 1e-7
    @test size(mesh.J) == (32, 16)
    @test maximum(abs.(vec(sum(mesh.J, dims = 2))[17:32] .- gasj(mesh)) ./ gasj(mesh)) < 1e-12

    # 2. anisotropic: conservative, non-negative, different from diffuse; GMRES matches a dense solve
    aniso = cavity([1000.0, 500.0, 500.0, 500.0]; phase = HenyeyGreenstein(0.8),
                   reflection = SpecularReflection(1.0))
    smooth!(aniso; verbose = false)
    solveEquilibrium!(aniso, aniso.F_smooth; verbose = false)
    @test abs(aniso.energy_error) < 1e-10
    @test minimum(aniso.J) > -1e-9 * maximum(aniso.J)
    @test maximum(abs.(gasT(aniso) .- T_diffuse)) > 1.0
    op, h0, _, _ = M._directional_setup(aniso, aniso.G_smooth)
    @test length(op.tables) == 5                                     # one HG table + four wall mirrors
    n = op.N * op.A
    Md = zeros(n, n); e = zeros(n); y = zeros(n)
    for c in 1:n
        fill!(e, 0.0); e[c] = 1.0
        mul!(y, op, e)
        Md[:, c] .= y
    end
    J_dense = Md \ vec(h0 .* op.S)
    @test norm(vec(aniso.J) .- J_dense) / norm(J_dense) < 1e-9

    # partial specularity sits between diffuse and mirror walls
    half = cavity([1000.0, 500.0, 500.0, 500.0]; reflection = SpecularReflection(0.5))
    smooth!(half; verbose = false)
    solveEquilibrium!(half, half.F_smooth; verbose = false)
    @test abs(half.energy_error) < 1e-10
    @test 0 < maximum(abs.(gasT(half) .- T_diffuse))

    # 3. isothermal enclosure: exact for any valid table once G is smoothed
    iso = cavity(fill(1000.0, 4); phase = HenyeyGreenstein(0.8), reflection = SpecularReflection(1.0))
    smooth!(iso; verbose = false)
    solveEquilibrium!(iso, iso.F_smooth; verbose = false)
    @test maximum(abs.(gasT(iso) .- 1000.0)) < 1e-6
    @test abs(iso.energy_error) < 1e-10
end

@testset "spectral directional solver" begin
    M = RayTraceHeatTransfer
    verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    dm = AngularBins(8, 2)
    K  = 3
    function spectral_cavity(T_walls; κ = [0.05, 1.0, 5.0], σs = [0.5, 0.8, 0.2], ε = [0.3, 0.5, 0.9],
                             phase = IsotropicScattering(), reflection = DiffuseReflection(),
                             model = PlanckBands([5e-7, 2e-6, 8e-6, 5e-5]))
        f = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), K, 1.0, 0.0)
        f.kappa_g = copy(κ); f.sigma_s_g = copy(σs); f.epsilon = [copy(ε) for _ in 1:4]
        f.T_in_w = copy(T_walls); f.q_in_w = zeros(4); f.T_in_g = -1.0; f.q_in_g = 0.0
        f.phase = phase; f.reflection = fill(reflection, 4)
        mesh = RayTracingDomain2D([f], [(3, 3)], verbose = false)
        mesh.spectral_mode = :spectral_variable          # also for bins with identical properties (test 3)
        mesh.spectral_model = model
        mesh.directional_model = dm
        mesh(200_000; method = :pathlength, seeds = 100, verbose = false)
        smooth!(mesh; verbose = false)
        return mesh
    end
    gasT(mesh) = [mesh.fine_mesh[ci][fi].T_g for ((ci, fi), _) in sort(collect(mesh.volume_mapping), by = last)]
    solve!(mesh) = solveEquilibrium!(mesh, mesh.F_smooth; max_iters = 5000, convergence_tol = 1e-12, verbose = false)

    # 1. default descriptors reproduce the spectral (diffuse) Woodbury solver on the same F_smooth
    spec = spectral_cavity([1000.0, 500.0, 500.0, 500.0])
    M.equilibriumSpectral2D!(spec, spec.F_smooth; max_iters = 5000, convergence_tol = 1e-12, verbose = false)
    T_diffuse = gasT(spec)
    solve!(spec)
    @test maximum(abs.(gasT(spec) .- T_diffuse)) < 1e-5
    @test maximum(abs.(spec.energy_error)) < 1e-9

    # 2. anisotropic, per-bin phase functions: conservative per bin, non-negative, different from diffuse
    aniso = spectral_cavity([1000.0, 500.0, 500.0, 500.0];
                            phase = [HenyeyGreenstein(0.2), HenyeyGreenstein(0.8), IsotropicScattering()],
                            reflection = SpecularReflection(0.7))
    solve!(aniso)
    @test maximum(abs.(aniso.energy_error)) < 1e-9
    @test length(aniso.J) == K && all(J -> size(J) == (21, 16), aniso.J)
    @test all(J -> minimum(J) > -1e-9 * maximum(J), aniso.J)
    @test maximum(abs.(gasT(aniso) .- T_diffuse)) > 0.5
    for k in 1:K                                          # Σ_b J = the radiosity written to the faces
        jv = [aniso.fine_mesh[ci][fi].j_g[k] for ((ci, fi), _) in sort(collect(aniso.volume_mapping), by = last)]
        @test maximum(abs.(vec(sum(aniso.J[k], dims = 2))[13:21] .- jv) ./ jv) < 1e-12
    end

    # 3. isothermal enclosure: exact for any valid tables
    iso = spectral_cavity(fill(1000.0, 4);
                          phase = [HenyeyGreenstein(0.2), HenyeyGreenstein(0.8), IsotropicScattering()],
                          reflection = SpecularReflection(1.0))
    solve!(iso)
    @test maximum(abs.(gasT(iso) .- 1000.0)) < 1e-5

    # 4. bins with identical properties and constant weights: the grey directional solution
    grey_like = spectral_cavity([1000.0, 500.0, 500.0, 500.0]; κ = fill(0.2, K), σs = fill(0.8, K), ε = fill(0.5, K),
                                phase = HenyeyGreenstein(0.8), reflection = SpecularReflection(1.0),
                                model = ConstantWeights([0.5, 0.3, 0.2]))
    solve!(grey_like)
    fg = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), 1, 0.2, 0.8)
    fg.epsilon = fill(0.5, 4); fg.T_in_w = [1000.0, 500.0, 500.0, 500.0]; fg.q_in_w = zeros(4)
    fg.T_in_g = -1.0; fg.q_in_g = 0.0
    fg.phase = HenyeyGreenstein(0.8); fg.reflection = fill(SpecularReflection(1.0), 4)
    grey = RayTracingDomain2D([fg], [(3, 3)], verbose = false)
    grey.directional_model = dm
    grey(200_000; method = :pathlength, seeds = 100, verbose = false)
    smooth!(grey; verbose = false)
    solveEquilibrium!(grey, grey.F_smooth; verbose = false)
    @test maximum(abs.(gasT(grey_like) .- gasT(grey))) < 1e-5
end

println("✓ Directional solvers tests complete")