# Directional models: angular bins, emission shares, redistribution tables.

println("\n" * "-"^60)
println("Testing directional models (bins, shares, tables)")
println("-"^60)

using Test
using RayTraceHeatTransfer
using Random, LinearAlgebra, SparseArrays
using GeometryBasics, StaticArrays

@testset "AngularBins" begin
    M = RayTraceHeatTransfer
    m = AngularBins(16, 4)
    @test M.n_angular_bins(m) == 64
    @test M.validate(m) === nothing
    @test_throws ErrorException M.validate(AngularBins(7, 2))
    @test_throws ErrorException M.validate(AngularBins(8, 0))
    for a in 1:64
        @test M.reversed_bin(m, a) != a
        @test M.reversed_bin(m, M.reversed_bin(m, a)) == a
    end
    # reversing a direction lands in the reversed bin (tracer convention: in-plane projection)
    rng = MersenneTwister(1)
    for _ in 1:2000
        θ = acos(1 - 2 * rand(rng)); ψ = 2π * rand(rng)
        d = (sin(θ) * cos(ψ), cos(θ))
        @test M.angular_bin(m, (-d[1], -d[2])) == M.reversed_bin(m, M.angular_bin(m, d))
    end
    @test M.angular_bin(m, (1.0, 0.0)) == 9                  # φ = 0, |μ_z| = 0: iφ = 9, iz = 1
    @test M.angular_bin(m, (0.0, 0.0)) > 48                  # along z: last polar ring
end

@testset "emission shares" begin
    M = RayTraceHeatTransfer
    m = AngularBins(16, 4)
    pv = M.emission_shares(m)
    @test all(pv .≈ 1 / 64) && sum(pv) ≈ 1
    for φn in (π / 2, 0.0, -π / 2, π, 0.3, -2.1)             # axis-aligned and oblique walls
        p = M.emission_shares(m, φn)
        @test sum(p) ≈ 1
        @test all(p .>= 0)
        # a bin and its reverse cannot both emit, except bins straddling the tangent
        n_both = count(a -> p[a] > 0 && p[M.reversed_bin(m, a)] > 0, 1:64)
        @test n_both <= 2 * 2 * m.n_polar
    end
    p = M.emission_shares(m, π / 2)                          # bottom wall, normal +y
    for iz in 1:4, iφ in 1:8                                 # φ ∈ [−π, 0): behind the wall
        @test p[(iz - 1) * 16 + iφ] == 0.0
    end
    # azimuthal marginal of a Lambertian wall: ½(sin φhi′ − sin φlo′) about the normal
    marg = [sum(p[(iz - 1) * 16 + iφ] for iz in 1:4) for iφ in 1:16]
    @test marg[13] ≈ 0.5 * (sin(π / 8) - sin(0)) atol = 1e-12   # bin [π/2, 5π/8): 0..π/8 off the normal
end

@testset "scattering tables" begin
    M = RayTraceHeatTransfer
    m = AngularBins(8, 2)
    A = 16
    @test M.redistribution_table(m, IsotropicScattering()) ≈ fill(1 / A, A, A)
    @test maximum(abs.(M.redistribution_table(m, HenyeyGreenstein(0.0)) .- 1 / A)) < 1e-13
    Φ = M.redistribution_table(m, HenyeyGreenstein(0.8))
    @test Φ ≈ Φ'
    @test maximum(abs.(sum(Φ, dims = 1) .- 1)) < 1e-12
    @test maximum(abs.(sum(Φ, dims = 2) .- 1)) < 1e-12
    for a in 1:A
        @test Φ[a, a] > Φ[a, M.reversed_bin(m, a)]           # forward-peaked
    end
    @test_throws ErrorException HenyeyGreenstein(1.0)
    @test M.redistribution_table(m, TabulatedScattering(Φ)) == Φ
    bad = copy(Φ); bad[1, :] .= 0; bad[1, 1] = 1.0           # rows still sum to 1, columns do not
    @test_throws ErrorException M.redistribution_table(m, TabulatedScattering(bad))
end

@testset "wall tables" begin
    M = RayTraceHeatTransfer
    m = AngularBins(16, 4)
    A = 64
    for φn in (π / 2, 0.0, -π / 2, π)                        # axis-aligned: exact permutation
        P = M.mirror_table(m, φn)
        p = M.emission_shares(m, φn)
        for b in 1:A
            p[M.reversed_bin(m, b)] > 0 || continue          # only rows that carry incident flux
            @test count(==(1.0), P[b, :]) == 1 && count(!=(0.0), P[b, :]) == 1
            bp = findfirst(==(1.0), P[b, :])
            @test p[bp] > 0                                  # mirrored into an emitting bin
            @test (bp - 1) ÷ 16 == (b - 1) ÷ 16              # |μ_z| ring kept
        end
    end
    for φn in (0.3, -2.1, 1.0), s in (1.0, 0.4, 0.0)         # oblique walls, partial specularity
        p = M.emission_shares(m, φn)
        Φ = M.redistribution_table(m, SpecularReflection(s), φn)   # checks itself on construction
        @test M.check_redistribution_table(Φ, p, m) === nothing
    end
    Φd = M.redistribution_table(m, DiffuseReflection(), 0.3)
    @test M.check_redistribution_table(Φd, M.emission_shares(m, 0.3), m) === nothing
    @test M.redistribution_table(m, SpecularReflection(0.0), 0.3) ≈ Φd
    @test_throws ErrorException SpecularReflection(1.5)
    @test_throws ErrorException M.redistribution_table(m, TabulatedReflection(fill(1 / A, A, A)), π / 2)
end

@testset "face descriptors: defaults and inheritance" begin

    quad_verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    tri_verts  = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(0.0, 1.0))

    # defaults reproduce the diffuse/isotropic behaviour, grey and spectral, quad and triangle
    for n in (1, 3)
        fq = PolyVolume2D{Float64}(quad_verts, SVector(true, true, true, true), n, 1.0, 0.0)
        ft = PolyVolume2D{Float64}(tri_verts, SVector(true, true, true), n, 1.0, 0.0)
        @test fq.phase === IsotropicScattering() && ft.phase === IsotropicScattering()
        @test length(fq.reflection) == 4 && length(ft.reflection) == 3
        @test all(r -> r === DiffuseReflection(), fq.reflection)
        @test all(r -> r === DiffuseReflection(), ft.reflection)
    end

    # which coarse wall does a fine wall lie on? 1 bottom, 2 right, 3 top, 4 left
    function coarse_wall(f, wi)
        p1 = f.vertices[wi]; p2 = f.vertices[mod1(wi + 1, length(f.vertices))]
        mx, my = (p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2
        return my < 1e-9 ? 1 : mx > 1 - 1e-9 ? 2 : my > 1 - 1e-9 ? 3 : 4
    end

    # grey: one descriptor per wall and per volume, inherited by every fine element
    face = PolyVolume2D{Float64}(quad_verts, SVector(true, true, true, true), 1, 0.2, 0.8)
    face.epsilon = fill(0.5, 4)
    face.T_in_w  = [1000.0, 500.0, 500.0, 500.0]
    face.q_in_w  = zeros(4)
    face.T_in_g  = -1.0
    face.q_in_g  = 0.0
    face.phase      = HenyeyGreenstein(0.8)
    face.reflection = [SpecularReflection(1.0), DiffuseReflection(), SpecularReflection(0.5), DiffuseReflection()]
    mesh = RayTracingDomain2D([face], [(3, 3)], verbose = false)
    @test all(f -> f.phase === HenyeyGreenstein(0.8), mesh.fine_mesh[1])
    @test length(mesh.surface_mapping) == 12
    for ((ci, fi, wi), _) in mesh.surface_mapping
        f = mesh.fine_mesh[ci][fi]
        @test f.reflection[wi] === face.reflection[coarse_wall(f, wi)]
    end

    # spectral: per-bin vectors are inherited as copies, scalar descriptors stay scalar
    K = 3
    sface = PolyVolume2D{Float64}(quad_verts, SVector(true, true, true, true), K, 1.0, 0.0)
    sface.kappa_g   = [0.1, 1.0, 5.0]
    sface.sigma_s_g = [0.5, 0.5, 0.5]
    sface.epsilon   = [fill(0.5, K) for _ in 1:4]
    sface.T_in_w    = [1000.0, 500.0, 500.0, 500.0]
    sface.q_in_w    = zeros(4)
    sface.T_in_g    = -1.0
    sface.q_in_g    = 0.0
    sface.phase      = [HenyeyGreenstein(0.1), HenyeyGreenstein(0.2), IsotropicScattering()]
    sface.reflection = [[SpecularReflection(1.0), SpecularReflection(0.5), DiffuseReflection()],
                        DiffuseReflection(), SpecularReflection(0.3), DiffuseReflection()]
    smesh = RayTracingDomain2D([sface], [(2, 2)], verbose = false)
    for f in smesh.fine_mesh[1]
        @test f.phase == sface.phase
        @test f.phase !== sface.phase                          # a copy, not an alias
    end
    for ((ci, fi, wi), _) in smesh.surface_mapping
        f = smesh.fine_mesh[ci][fi]
        cw = coarse_wall(f, wi)
        if cw == 1
            @test f.reflection[wi] == sface.reflection[1]
            @test f.reflection[wi] !== sface.reflection[1]
        else
            @test f.reflection[wi] === sface.reflection[cw]
        end
    end
end

@testset "domain: directional model field, G storage, show" begin
    verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    face  = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), 1, 0.2, 0.8)
    face.epsilon = fill(0.5, 4)
    face.T_in_w  = [1000.0, 500.0, 500.0, 500.0]
    face.q_in_w  = zeros(4)
    face.T_in_g  = -1.0
    face.q_in_g  = 0.0
    mesh = RayTracingDomain2D([face], [(3, 3)], verbose = false)

    # defaults: a domain without a directional model is exactly the old domain
    @test mesh.directional_model === nothing
    @test mesh.G_raw === nothing && mesh.G_smooth === nothing
    txt = sprint(show, MIME"text/plain"(), mesh)
    @test !occursin("direction", txt)
    @test occursin("method = :exchange", txt)

    # with a model: reported, and the hint points at the pathlength tracer
    mesh.directional_model = AngularBins(16, 4)
    txt = sprint(show, MIME"text/plain"(), mesh)
    @test occursin("direction  AngularBins: 16 azimuthal × 4 out-of-plane = 64 direction bins", txt)
    @test occursin("method = :pathlength", txt)
    @test sprint(show, mesh) isa String                       # one-line form still prints

    # after an ordinary trace the G lines appear and say what is missing
    mesh(200_000; method = :pathlength, verbose = false)
    txt = sprint(show, MIME"text/plain"(), mesh)
    @test occursin("G_raw     64 × (", txt)
    @test occursin("G_smooth  not computed", txt)
    @test length(mesh.G_raw) == 64

    # the storage accepts both layouts
    N = length(mesh.surface_mapping) + length(mesh.volume_mapping)
    mesh.G_raw = [spzeros(N, N) for _ in 1:64]
    @test occursin("G_raw     64 × (", sprint(show, MIME"text/plain"(), mesh))
    mesh.G_raw = [[spzeros(N, N) for _ in 1:64] for _ in 1:3]
    @test occursin("G_raw     3 × (64 × (", sprint(show, MIME"text/plain"(), mesh))
    mesh.G_raw = nothing
end

@testset "angular deposition" begin
    M = RayTraceHeatTransfer
    verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    function cavity(κ, σs)
        n = length(κ)
        f = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), n, 1.0, 0.0)
        if n == 1
            f.kappa_g = κ[1]; f.sigma_s_g = σs[1]; f.epsilon = fill(0.5, 4)
        else
            f.kappa_g = copy(κ); f.sigma_s_g = copy(σs); f.epsilon = [fill(0.5, n) for _ in 1:4]
        end
        f.T_in_w = [1000.0, 500.0, 500.0, 500.0]; f.q_in_w = zeros(4)
        f.T_in_g = -1.0; f.q_in_g = 0.0
        return RayTracingDomain2D([f], [(4, 4)], verbose = false)
    end
    n_rays = 400_000
    dm = AngularBins(8, 2)
    A  = 16

    # grey: same seeds with and without a model give the same F_raw
    plain = cavity([0.2], [0.8])
    plain(n_rays; method = :pathlength, seeds = 100, verbose = false)
    @test plain.G_raw === nothing
    mesh = cavity([0.2], [0.8]); mesh.directional_model = dm
    mesh(n_rays; method = :pathlength, seeds = 100, verbose = false)
    G = mesh.G_raw
    N = size(mesh.F_raw, 1)
    ns = length(mesh.surface_mapping)
    @test G isa Vector && length(G) == A && all(g -> size(g) == (N, N), G)
    @test all(g -> all(>=(0), nonzeros(g)), G)
    @test maximum(abs.(sum(G) - mesh.F_raw)) < 1e-15                    # the invariant
    @test norm(mesh.F_raw - plain.F_raw) / norm(plain.F_raw) < 1e-12    # same F as the diffuse path
    @test maximum(abs.(sum(mesh.F_raw, dims = 2) .- 1)) < 1e-12
    @test mesh.G_smooth === nothing

    # emission shares: volumes isotropic within noise, walls empty behind themselves
    share = hcat([vec(sum(g, dims = 2)) for g in G]...)                 # N × A
    n_e = n_rays ÷ N
    σ = sqrt((1 / A) * (1 - 1 / A) / n_e)
    @test maximum(abs.(share[ns+1:end, :] .- 1 / A)) < 6σ
    for s in 1:ns
        @test count(==(0.0), share[s, :]) == A ÷ 2                      # axis-aligned walls, Aφ multiple of 4
    end

    # chunked trace: paths not kept, and the same F as a plain chunked trace (same seeds, same
    # chunking). A chunked trace visits the emitters once per chunk, so it consumes the per-thread
    # random streams in a different order than a single-chunk trace: those two agree only statistically.
    chunked = cavity([0.2], [0.8]); chunked.directional_model = dm
    chunked(n_rays; method = :pathlength, seeds = 100, chunk_rays = n_rays ÷ 4, verbose = false)
    chunked_plain = cavity([0.2], [0.8])
    chunked_plain(n_rays; method = :pathlength, seeds = 100, chunk_rays = n_rays ÷ 4, verbose = false)
    @test chunked.path_store === nothing
    @test length(chunked.G_raw) == A
    @test maximum(abs.(sum(chunked.G_raw) - chunked.F_raw)) < 1e-15
    @test norm(chunked.F_raw - chunked_plain.F_raw) / norm(chunked_plain.F_raw) < 1e-12
    @test maximum(maximum(abs.(chunked.G_raw[a] - G[a])) for a in 1:A) < 0.03   # statistical agreement only

    # re-binning from the kept store: new κ, new G, invariant kept; dropping the model clears G
    for f in mesh.fine_mesh[1]
        f.kappa_g = 2.0
    end
    exchangeFactors!(mesh; verbose = false)
    @test maximum(maximum(abs.(mesh.G_raw[a] - G[a])) for a in 1:A) > 1e-3
    @test maximum(abs.(sum(mesh.G_raw) - mesh.F_raw)) < 1e-15
    mesh.directional_model = nothing
    exchangeFactors!(mesh; verbose = false)
    @test mesh.G_raw === nothing

    # spectral layout: G_raw[k][a], invariant per spectral bin
    spec = cavity([0.05, 1.0, 5.0], [0.5, 0.5, 0.5]); spec.directional_model = dm
    spec(n_rays; method = :pathlength, seeds = 100, verbose = false)
    @test length(spec.G_raw) == 3 && all(g -> length(g) == A, spec.G_raw)
    for k in 1:3
        @test maximum(abs.(sum(spec.G_raw[k]) - spec.F_raw[k])) < 1e-15
        @test maximum(abs.(sum(spec.F_raw[k], dims = 2) .- 1)) < 1e-12
    end
end

@testset "smoothing of angular exchange factors" begin
    M = RayTraceHeatTransfer
    verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    dm = AngularBins(8, 2)
    A  = 16

    function check_smoothed(G, F, w, p)
        N = size(F, 1)
        @test all(g -> all(>=(0), nonzeros(g)), G)
        share = hcat([vec(sum(g, dims = 2)) for g in G]...)
        @test maximum(abs.(share .- p)) < 1e-12                          # per-bin emission shares
        recip = maximum(maximum(abs.(Diagonal(w) * G[a] - (Diagonal(w) * G[M.reversed_bin(dm, a)])'))
                        for a in 1:A)
        @test recip / maximum(w) < 1e-12                                 # reversed-bin reciprocity
        @test maximum(abs.(sum(G) - F)) < 1e-14                          # F_smooth = Σₐ G_smooth
        @test maximum(abs.(sum(F, dims = 2) .- 1)) < 1e-12
        WF = Diagonal(w) * Matrix(F)
        @test norm(WF - WF') / norm(WF) < 1e-12                          # ordinary reciprocity
    end

    # grey
    f = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), 1, 0.2, 0.8)
    f.epsilon = fill(0.5, 4); f.T_in_w = [1000.0, 500.0, 500.0, 500.0]; f.q_in_w = zeros(4)
    f.T_in_g = -1.0; f.q_in_g = 0.0
    mesh = RayTracingDomain2D([f], [(4, 4)], verbose = false)
    mesh.directional_model = dm
    mesh(400_000; method = :pathlength, seeds = 100, verbose = false)
    stats = smooth!(mesh; verbose = false)
    @test all(stats.converged)
    @test stats.k_dykstra == [0] && stats.k_pcg_tot == [0]
    @test stats.k_ap[1] > 0 && stats.delta_smooth[1] < 1e-10
    @test length(mesh.G_smooth) == A
    p = M.emission_share_matrix(mesh)
    check_smoothed(mesh.G_smooth, mesh.F_smooth, M.get_w(mesh), p)
    for a in 1:A                                                         # zero pattern kept
        @test all(mesh.G_raw[a][i, j] > 0 for (i, j, _) in zip(findnz(mesh.G_smooth[a])...))
    end
    @test occursin("G_smooth  16 × (", sprint(show, MIME"text/plain"(), mesh))

    # a k_ap that is too small is reported, not hidden
    stats_short = @test_logs (:warn,) match_mode = :any smooth!(mesh; k_ap = 5, verbose = false)
    @test !all(stats_short.converged)

    # a new deposition invalidates the smoothed G
    smooth!(mesh; verbose = false)
    exchangeFactors!(mesh; verbose = false)
    @test mesh.G_smooth === nothing

    # spectral layout
    K = 3
    fs = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), K, 1.0, 0.0)
    fs.kappa_g = [0.05, 1.0, 5.0]; fs.sigma_s_g = [0.5, 0.5, 0.5]
    fs.epsilon = [fill(0.5, K) for _ in 1:4]; fs.T_in_w = [1000.0, 500.0, 500.0, 500.0]
    fs.q_in_w = zeros(4); fs.T_in_g = -1.0; fs.q_in_g = 0.0
    spec = RayTracingDomain2D([fs], [(4, 4)], verbose = false)
    spec.directional_model = dm
    spec(400_000; method = :pathlength, seeds = 100, verbose = false)
    sstats = smooth!(spec; verbose = false)
    @test length(sstats.converged) == K && all(sstats.converged)
    ps = M.emission_share_matrix(spec)
    for k in 1:K
        check_smoothed(spec.G_smooth[k], spec.F_smooth[k], M.get_w(spec; spectral_bin = k), ps)
    end
end

@testset "dense spectral-directional blocks" begin
    M = RayTraceHeatTransfer
    rng = MersenneTwister(7)
    verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    K = 3
    f = PolyVolume2D{Float64}(verts, SVector(true, true, true, true), K, 1.0, 0.0)
    f.kappa_g = [0.05, 1.0, 5.0]; f.sigma_s_g = [0.5, 0.8, 0.2]
    f.epsilon = [[0.3, 0.5, 0.9] for _ in 1:4]
    f.T_in_w = [1000.0, 500.0, -1.0, 500.0]; f.q_in_w = zeros(4)     # one flux-specified wall
    f.T_in_g = -1.0; f.q_in_g = 0.0
    f.phase = [HenyeyGreenstein(0.2), HenyeyGreenstein(0.8), IsotropicScattering()]   # per spectral bin
    f.reflection = fill(SpecularReflection(0.7), 4)                                   # all bins
    spec = RayTracingDomain2D([f], [(3, 3)], verbose = false)
    spec.directional_model = AngularBins(8, 2)
    spec(200_000; method = :pathlength, seeds = 100, verbose = false)
    smooth!(spec; verbose = false)
    b_pkg = M.get_b(spec)

    for k in 1:K
        op, _, _, b, known_q = M._directional_setup(spec, spec.G_smooth[k]; spectral_bin = k)
        N, A = op.N, op.A
        n = N * A
        @test b ≈ b_pkg[:, k]                                         # same albedo as the spectral solver
        @test count(known_q) == 3 + 9                                 # top wall elements + all volumes
        coeff = ifelse.(known_q, 1.0, b)
        D, Mb = M._directional_dense_blocks(op, coeff)
        @test size(D) == (n, n) && size(Mb) == (N, n)

        # 𝐃 equals the matrix-free operator with the re-emission term switched off
        op0 = M.DirectionalSystemOp(op.N, op.A, op.G, op.S, op.invS, op.refl, zeros(N), op.diff_w,
                                    op.tab_w, op.tab_idx, op.tables, zeros(N, A), zeros(N, A), zeros(n))
        x = rand(rng, n); y = zeros(n)
        mul!(y, op0, x)
        @test norm(D * x - y) / norm(y) < 1e-13

        # 𝐌: total outgoing − coeff · total incident
        M._incident!(op0, reshape(x, N, A))
        g_tot = vec(sum(op0.g, dims = 2))
        @test norm(Mb * x - (M._sum_bins(op, x) .- coeff .* g_tot)) / norm(x) < 1e-13

        # P and R, and the collapse onto the diffuse blocks of the spectral solver
        P = zeros(n, N); R = zeros(N, n)
        for a in 1:A, i in 1:N
            P[(a - 1) * N + i, i] = op.S[i, a]
            R[i, (a - 1) * N + i] = 1.0
        end
        e = rand(rng, N)
        @test M._spread_emission(op, e) ≈ P * e
        @test M._sum_bins(op, x) ≈ R * x
        F = Matrix(spec.F_smooth[k])
        @test maximum(abs.(R * D * P - (I - Diagonal(b) * F'))) < 1e-12
        @test maximum(abs.(Mb * P - (I - Diagonal(coeff) * F'))) < 1e-12
        @test maximum(abs.(Mb - R * D)[.!known_q, :]) < 1e-12          # border = summed block rows where coeff = ρ
    end
end

println("✓ Directional models tests complete")