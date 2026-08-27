
"""
Regression test (3D): grey vs spectral consistency with nonuniform b.

Pins the 3D sibling of the receiver-vs-sender b-indexing bug fixed July
2026: equilibriumSurfacesGrey3D! Step 8 computed absorbed/reflected energy
with B_k (sender's reflectivity) instead of B_i (receiver's). Wrong r =>
wrong e = j - r => wrong recovered temperatures on radiative-equilibrium
faces, whenever epsilon is NONUNIFORM across faces. Every prior 3D test
used uniform epsilon and was blind to it.

Construction: cube with mixed per-face emissivity, two prescribed-T faces,
four radiative-equilibrium faces. The spectral twin replicates each face's
epsilon uniformly over 3 bins, so the spectral solve must reproduce the
grey solve exactly (bin-uniform properties => weighted and plain emission
fractions coincide; the spectral problem is the grey problem partitioned).
3D view factors are deterministic, so both domains produce the same F and
the comparison carries no Monte Carlo allowance: solver-precision
tolerances. A regression fails by tens of percent, not marginally.
"""

println("\n" * "-"^60)
println("Testing 3D Grey vs Spectral Consistency (Nonuniform b)")
println("-"^60)

using RayTraceHeatTransfer
using Test
using LinearAlgebra
using StatsBase
using StaticArrays
using GeometryBasics
using SparseArrays
using ConvolutionInterpolations

GREY3D_TEMP_TOL = 1e-6   # K; deterministic F, solver-precision comparison
GREY3D_J_RTOL   = 1e-9

@testset "3D Shared Grey vs Spectral (mixed epsilon)" begin
    # ---- unit cube, standard connectivity from the existing 3D tests ---------
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

    Ndim   = 4
    n_bins = 3

    # Two prescribed faces, four radiative-equilibrium faces
    T_in_w = [1000.0, 300.0, -1.0, -1.0, -1.0, -1.0]
    q_in_w = zeros(6)

    # MIXED emissivity: the configuration the B_k bug corrupts.
    # Faces 1-3 reflective (b = 0.7), faces 4-6 nearly black (b = 0.1).
    eps_faces = [0.3, 0.3, 0.3, 0.9, 0.9, 0.9]

    # ---- grey domain ---------------------------------------------------------
    domain_grey = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w,
                                     copy(eps_faces))
    domain_grey()
    smooth!(domain_grey)
    solveEquilibrium!(domain_grey, domain_grey.F_smooth)

    # ---- spectral twin: same epsilon replicated over bins --------------------
    eps_spectral = [fill(eps_faces[i], n_bins) for i in 1:6]
    domain_spec = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w,
                                     eps_spectral)
    domain_spec()
    smooth!(domain_spec)
    domain_spec.wavelength_band_limits = [1.0e-6, 3.0e-6, 8.0e-6, 1.0e-3]
    solveEquilibrium!(domain_spec, domain_spec.F_smooth)

    # ---- collect per-subface T and j -----------------------------------------
    T_grey = Float64[]; j_grey = Float64[]
    for superface in domain_grey.facesMesh
        for subface in superface.subFaces
            push!(T_grey, subface.T_w)
            push!(j_grey, subface.j_w isa AbstractVector ? sum(subface.j_w) :
                                                           subface.j_w)
        end
    end

    T_spec = Float64[]; j_spec = Float64[]
    for superface in domain_spec.facesMesh
        for subface in superface.subFaces
            push!(T_spec, subface.T_w)
            jw = subface.j_w
            @test jw isa AbstractVector       # spectral shape asserted
            @test length(jw) == n_bins
            push!(j_spec, sum(jw))
        end
    end

    @test length(T_grey) == length(T_spec)

    # ---- temperatures: must agree to solver precision ------------------------
    # The B_k bug produced wrong r on every subface receiving from faces of
    # differing epsilon, corrupting recovered T on the equilibrium faces.
    @test maximum(abs.(T_spec .- T_grey)) < GREY3D_TEMP_TOL

    # ---- radiosity: bin-summed spectral must equal grey ----------------------
    @test all(isapprox.(j_spec, j_grey; rtol=GREY3D_J_RTOL, atol=1e-8))

    # ---- direct pin on the post-processing split -----------------------------
    # Per subface: e = j - r must satisfy e = eps*sigma*A*T^4 with the
    # RECEIVER's eps (grey path), equivalently the recovered T must be
    # consistent with the receiver-indexed split. Cross-checking T^4 against
    # the bin-summed spectral e (which uses D = I - Diag(b)F' correctly)
    # pins the grey Step 8 specifically.
    e_spec = Float64[]
    for superface in domain_spec.facesMesh
        for subface in superface.subFaces
            push!(e_spec, sum(subface.e_w))
        end
    end
    e_grey = Float64[]
    for superface in domain_grey.facesMesh
        for subface in superface.subFaces
            ew = subface.e_w
            push!(e_grey, ew isa AbstractVector ? sum(ew) : ew)
        end
    end
    mask = e_grey .> 1e-8
    @test all(isapprox.(e_spec[mask], e_grey[mask]; rtol=1e-9))
end

println("✓ 3D Grey vs Spectral consistency tests complete")