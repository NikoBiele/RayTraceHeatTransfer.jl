# test_spectral_bin_refinement.jl
#
# Nested-bin refinement invariance: splitting every spectral bin in two,
# with children inheriting the parent bin's properties, leaves the physics
# EXACTLY unchanged — Planck fractions telescope (f_parent = f_child1 +
# f_child2 identically), so the solved temperature field must agree to
# solver precision, and per-bin radiosities must satisfy
# j_parent = j_child1 + j_child2.
#
# This is an exact physics oracle for the entire spectral machinery with
# genuinely bin-varying properties: emission-fraction computation and
# weighting, per-bin system assembly, bin bookkeeping/indexing, and the
# solution writeback. An off-by-one in bin indexing or a wrong weight
# normalization fails this test loudly while passing grey-vs-spectral
# comparisons (which require spectrally uniform properties to have a grey
# reference).
#
# CRITICAL DESIGN POINT: both solves must run on the SAME exchange matrices.
# In 2D the coarse trace's F_smooth (one matrix per distinct extinction) is
# duplicated onto the child bins — the fine mesh is NOT re-traced, so Monte
# Carlo noise cancels exactly and the ray count can stay modest. In 3D the
# view factors are deterministic and band-independent, so both domains
# compute identical F; this is asserted as a precondition.
#
# Tolerances: solves are run to convergence_tol = 1e-14; the comparison
# floor is solver arithmetic (~1e-11 K, post-v0.12.3 unridged), tested at
# 1e-9 K / rtol 1e-8 with the usual ~100x margin.

using Test
using RayTraceHeatTransfer
using GeometryBasics, StaticArrays

REFINE_T_TOL  = 1e-9   # K
REFINE_J_RTOL = 1e-8

"Refine band edges by inserting the geometric midpoint of every bin."
function refine_edges(edges::Vector{Float64})
    fine = Float64[edges[1]]
    for i in 1:length(edges)-1
        push!(fine, sqrt(edges[i] * edges[i+1]))
        push!(fine, edges[i+1])
    end
    return fine
end

"Duplicate per-bin values onto child bins: [a, b, c] -> [a, a, b, b, c, c]."
duplicate_bins(v::Vector) = reduce(vcat, [[x, x] for x in v])

@testset "2D nested-bin refinement invariance" begin
    K_coarse = 4
    K_fine   = 2 * K_coarse
    Ndim     = 4
    N_rays   = 400_000        # modest on purpose: shared F, oracle is exact

    edges_coarse = 10 .^ range(log10(1e-7), log10(1e-3), length = K_coarse + 1)
    edges_fine   = refine_edges(collect(edges_coarse))

    # Genuinely bin-varying properties, with reflection AND scattering (b > 0)
    kappa_c  = [0.4, 0.8, 1.2, 1.6]
    sigma_c  = [0.2, 0.1, 0.3, 0.15]
    eps_c    = [0.9, 0.7, 0.8, 0.6]

    function build_2d(K, edges, kappa_v, sigma_v, eps_v)
        verts = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0),
                        Point2(1.0, 1.0), Point2(0.0, 1.0))
        solid = SVector(true, true, true, true)
        face  = PolyVolume2D{Float64}(verts, solid, K, kappa_v[1], sigma_v[1])
        face.kappa_g   = copy(kappa_v)
        face.sigma_s_g = copy(sigma_v)
        face.epsilon   = [copy(eps_v) for _ in 1:4]
        face.T_in_w    = [1000.0, 0.0, 0.0, 0.0]     # anchored: unique solution
        face.q_in_w    = zeros(4)
        face.T_in_g    = -1.0
        face.q_in_g    = 0.0
        mesh = RayTracingDomain2D([face], [(Ndim, Ndim)])
        mesh.wavelength_band_limits = edges
        return mesh
    end

    # ---- coarse: trace, smooth, solve ---------------------------------------
    mesh_c = build_2d(K_coarse, collect(edges_coarse), kappa_c, sigma_c, eps_c)
    mesh_c(N_rays; method = :exchange)
    smooth!(mesh_c; verbose = false)
    solveEquilibrium!(mesh_c, mesh_c.F_smooth;
                      max_iters = 10_000, convergence_tol = 1e-14, verbose = false)
    T_c = [ff.T_g for ff in mesh_c.fine_mesh[1]]

    # ---- fine: SAME F, duplicated per child bin — no re-trace ---------------
    F_fine = AbstractMatrix[mesh_c.F_smooth[cld(k, 2)] for k in 1:K_fine]
    mesh_f = build_2d(K_fine, edges_fine,
                      duplicate_bins(kappa_c), duplicate_bins(sigma_c),
                      duplicate_bins(eps_c))
    mesh_f.F_smooth = F_fine       # intended bypass of the trace; F is shared
    solveEquilibrium!(mesh_f, F_fine;
                      max_iters = 10_000, convergence_tol = 1e-14, verbose = false)
    T_f = [ff.T_g for ff in mesh_f.fine_mesh[1]]

    # ---- exact-physics oracle: temperatures unchanged by refinement ---------
    @test length(T_c) == length(T_f)
    @test maximum(abs.(T_c .- T_f)) < REFINE_T_TOL
end

@testset "3D nested-bin refinement invariance" begin
    # Unit cube, standard connectivity from the existing 3D tests
    points = [
        0.0 0.0 0.0;
        0.0 0.0 1.0;
        0.0 1.0 0.0;
        0.0 1.0 1.0;
        1.0 0.0 0.0;
        1.0 0.0 1.0;
        1.0 1.0 0.0;
        1.0 1.0 1.0
    ]
    faces = [
        1 2 4 3;
        5 6 8 7;
        1 5 7 3;
        2 6 8 4;
        3 4 8 7;
        1 2 6 5
    ]

    Ndim     = 4
    nb_c     = 3
    nb_f     = 2 * nb_c

    T_in_w = [1000.0, 300.0, -1.0, -1.0, -1.0, -1.0]
    q_in_w = zeros(6)

    edges_c = [1.0e-6, 3.0e-6, 8.0e-6, 1.0e-3]
    edges_f = refine_edges(edges_c)

    # Strongly bin-varying emissivity per face (b > 0 everywhere, no grey twin)
    eps_c_faces = [[0.3, 0.7, 0.5],
                   [0.9, 0.4, 0.6],
                   [0.5, 0.8, 0.3],
                   [0.7, 0.2, 0.9],
                   [0.4, 0.6, 0.8],
                   [0.8, 0.5, 0.4]]

    # ---- coarse -------------------------------------------------------------
    dom_c = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w,
                               [copy(e) for e in eps_c_faces])
    dom_c()
    smooth!(dom_c; verbose = false)
    dom_c.wavelength_band_limits = edges_c
    solveEquilibrium!(dom_c, dom_c.F_smooth;
                      max_iters = 10_000, convergence_tol = 1e-14, verbose = false)

    # ---- fine: children inherit the parent bin's epsilon --------------------
    dom_f = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w,
                               [duplicate_bins(e) for e in eps_c_faces])
    dom_f()
    smooth!(dom_f; verbose = false)
    dom_f.wavelength_band_limits = edges_f

    # Precondition: deterministic band-independent view factors => identical F
    @test maximum(abs.(Matrix(dom_f.F_smooth) .- Matrix(dom_c.F_smooth))) < 1e-14

    solveEquilibrium!(dom_f, dom_f.F_smooth;
                      max_iters = 10_000, convergence_tol = 1e-14, verbose = false)

    # ---- collect and compare ------------------------------------------------
    T_c = Float64[]; J_c = Vector{Float64}[]
    for sf in dom_c.facesMesh, sub in sf.subFaces
        push!(T_c, sub.T_w); push!(J_c, copy(sub.j_w))
    end
    T_f = Float64[]; J_f = Vector{Float64}[]
    for sf in dom_f.facesMesh, sub in sf.subFaces
        push!(T_f, sub.T_w); push!(J_f, copy(sub.j_w))
    end

    @test length(T_c) == length(T_f)
    @test maximum(abs.(T_c .- T_f)) < REFINE_T_TOL

    # Per-bin radiosity telescoping: j_parent = j_child1 + j_child2
    for i in eachindex(J_c)
        @test length(J_c[i]) == nb_c && length(J_f[i]) == nb_f
        for k in 1:nb_c
            jc = J_c[i][k]
            jf = J_f[i][2k-1] + J_f[i][2k]
            @test isapprox(jc, jf; rtol = REFINE_J_RTOL,
                                   atol = REFINE_J_RTOL * max(1.0, abs(jc)))
        end
    end
end