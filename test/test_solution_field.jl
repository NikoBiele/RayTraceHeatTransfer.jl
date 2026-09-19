# The global solution field J: present on every solver path, in global element
# order, with the layout of the exchange factors it pairs with.

println("\n" * "-"^60)
println("Testing the global solution field J")
println("-"^60)

using Test
using RayTraceHeatTransfer
using GeometryBasics, StaticArrays, LinearAlgebra, SparseArrays

@testset "J on the 2D paths" begin
    verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    # element vector in global order from a wall accessor and a volume accessor
    function elementvec(mesh, fw, fg)
        ns = length(mesh.surface_mapping)
        v  = zeros(ns + length(mesh.volume_mapping))
        for ((ci, fi, wi), s) in mesh.surface_mapping
            v[s] = fw(mesh.fine_mesh[ci][fi], wi)
        end
        for ((ci, fi), k) in mesh.volume_mapping
            v[ns + k] = fg(mesh.fine_mesh[ci][fi])
        end
        return v
    end

    # grey: J is the radiosity vector, and Fᵀ J is the incident power written to the faces
    f = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), 1, 0.2, 0.8)
    f.epsilon = fill(0.5, 4); f.T_in_w = [1000.0, 500.0, 500.0, 500.0]; f.q_in_w = zeros(4)
    f.T_in_g = -1.0; f.q_in_g = 0.0
    grey = RayTracingDomain2D([f], [(4, 4)], verbose = false)
    @test grey.J === nothing
    grey(400_000; method = :pathlength, seeds = 100, verbose = false)
    smooth!(grey; verbose = false)
    solveEquilibrium!(grey, grey.F_smooth; verbose = false)
    @test grey.J isa Vector{Float64} && length(grey.J) == 32
    @test grey.J ≈ elementvec(grey, (ff, wi) -> ff.j_w[wi], ff -> ff.j_g)
    @test grey.F_smooth' * grey.J ≈ elementvec(grey, (ff, wi) -> ff.g_w[wi], ff -> ff.g_g)

    # grey directional: J[i, b]; its sum over b is the radiosity, and the incident power follows from
    # G with the share normalisation: g[:, b] = G[b]ᵀ (J[:, b] ./ S[:, b]), S = row sums of G[b]
    fd = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), 1, 0.2, 0.8)
    fd.epsilon = fill(0.5, 4); fd.T_in_w = [1000.0, 500.0, 500.0, 500.0]; fd.q_in_w = zeros(4)
    fd.T_in_g = -1.0; fd.q_in_g = 0.0
    fd.phase = HenyeyGreenstein(0.8); fd.reflection = fill(SpecularReflection(1.0), 4)
    dir = RayTracingDomain2D([fd], [(4, 4)], verbose = false)
    dir.directional_model = AngularBins(8, 2)
    dir(400_000; method = :pathlength, seeds = 100, verbose = false)
    smooth!(dir; verbose = false)
    solveEquilibrium!(dir, dir.F_smooth; verbose = false)
    @test dir.J isa Matrix{Float64} && size(dir.J) == (32, 16)
    @test vec(sum(dir.J, dims = 2)) ≈ elementvec(dir, (ff, wi) -> ff.j_w[wi], ff -> ff.j_g)
    g = zeros(32)
    for b in 1:16
        S = vec(sum(dir.G_smooth[b], dims = 2))
        g .+= dir.G_smooth[b]' * (dir.J[:, b] .* [s > 0 ? 1 / s : 0.0 for s in S])
    end
    @test g ≈ elementvec(dir, (ff, wi) -> ff.g_w[wi], ff -> ff.g_g)

    # spectral: J[k][i]
    K = 3
    fs = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), K, 1.0, 0.0)
    fs.kappa_g = [0.05, 1.0, 5.0]; fs.sigma_s_g = [0.5, 0.8, 0.2]
    fs.epsilon = [[0.3, 0.5, 0.9] for _ in 1:4]; fs.T_in_w = [1000.0, 500.0, 500.0, 500.0]
    fs.q_in_w = zeros(4); fs.T_in_g = -1.0; fs.q_in_g = 0.0
    spec = RayTracingDomain2D([fs], [(3, 3)], verbose = false)
    spec.spectral_model = PlanckBands([5e-7, 2e-6, 8e-6, 5e-5])
    spec(200_000; method = :pathlength, seeds = 100, verbose = false)
    smooth!(spec; verbose = false)
    solveEquilibrium!(spec, spec.F_smooth; max_iters = 2000, convergence_tol = 1e-12, verbose = false)
    @test spec.J isa Vector{Vector{Float64}} && length(spec.J) == K
    for k in 1:K
        @test spec.J[k] ≈ elementvec(spec, (ff, wi) -> ff.j_w[wi][k], ff -> ff.j_g[k])
        @test spec.F_smooth[k]' * spec.J[k] ≈ elementvec(spec, (ff, wi) -> ff.g_w[wi][k], ff -> ff.g_g[k])
    end

    # spectral directional: J[k][i, b]
    spec.directional_model = AngularBins(8, 2)
    exchangeFactors!(spec; verbose = false)
    smooth!(spec; verbose = false)
    solveEquilibrium!(spec, spec.F_smooth; max_iters = 5000, convergence_tol = 1e-12, verbose = false)
    @test spec.J isa Vector{Matrix{Float64}} && length(spec.J) == K
    for k in 1:K
        @test vec(sum(spec.J[k], dims = 2)) ≈ elementvec(spec, (ff, wi) -> ff.j_w[wi][k], ff -> ff.j_g[k])
    end
end

@testset "J on the 3D surface paths" begin
    points = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0; 0.0 1.0 1.0;
              1.0 0.0 0.0; 1.0 0.0 1.0; 1.0 1.0 0.0; 1.0 1.0 1.0]
    faces  = [1 2 4 3; 5 6 8 7; 1 5 7 3; 2 6 8 4; 3 4 8 7; 1 2 6 5]
    T_in_w = [1000.0, 500.0, 500.0, 500.0, 500.0, 500.0]
    q_in_w = zeros(6)

    grey = ViewFactorDomain3D(points, faces, 2, q_in_w, T_in_w, fill(0.7, 6))
    @test grey.J === nothing
    grey(; verbose = false)
    smooth!(grey, verbose = false)
    solveEquilibrium!(grey, grey.F_smooth, verbose = false)
    @test grey.J isa Vector{Float64} && length(grey.J) == 24
    @test grey.J ≈ [sf.j_w for sf in elementFaces(grey)]

    # two spectral emissivity layouts, so that whichever 3D spectral solver each is routed to is covered
    for eps in ([[0.3, 0.6, 0.9], [0.5, 0.6, 0.7], [0.5, 0.6, 0.7], [0.5, 0.6, 0.7], [0.5, 0.6, 0.7], [0.5, 0.6, 0.7]],
                [[0.5, 0.6, 0.7] for _ in 1:6])
        spec = ViewFactorDomain3D(points, faces, 2, q_in_w, T_in_w, eps)
        spec.spectral_model = PlanckBands([5e-7, 2e-6, 8e-6, 5e-5])
        spec(; verbose = false)
        smooth!(spec, verbose = false)
        solveEquilibrium!(spec, spec.F_smooth, verbose = false)
        @test spec.J isa Vector{Vector{Float64}} && length(spec.J) == 3
        for k in 1:3
            @test spec.J[k] ≈ [sf.j_w[k] for sf in elementFaces(spec)]
        end
    end
end

println("✓ Global solution field J tests complete")