# An exactly transparent spectral bin (κ = 0) in a participating domain.
#
# Guards three fixes found together:
#   - get_w no longer floors volume reciprocity weights at 1e-6; zero-weight
#     elements are excluded from the smoothing projection (DkAP_active)
#   - get_b returns zero albedo for β = 0 instead of 0/0
#   - the tracer handles β = 0 (infinite free path)
#
# Isothermal enclosure: T_g must equal T_wall to solver precision for black
# and reflecting walls, once with PlanckBands and once with ConstantWeights
# (the latter's only remaining coverage). Also checks the smoothed F of the
# transparent bin has exactly zero volume columns.

println("\n" * "-"^60)
println("Testing an exactly transparent spectral bin")
println("-"^60)

using Test
using RayTraceHeatTransfer
using GeometryBasics, StaticArrays

function _transparent_face(κ_bins, ε_wall, T_wall)
    n_bins = length(κ_bins)
    vertices = SVector(Point2(0.0, 0.0), Point2(1.0, 0.0), Point2(1.0, 1.0), Point2(0.0, 1.0))
    face = PolyVolume2D{Float64}(vertices, SVector(true, true, true, true), n_bins, 1.0, 0.0)
    face.kappa_g   = copy(κ_bins)
    face.sigma_s_g = zeros(n_bins)
    face.epsilon   = [fill(ε_wall, n_bins) for _ in 1:4]
    face.T_in_w    = fill(T_wall, 4)
    face.q_in_w    = zeros(4)
    face.T_in_g    = -1.0
    face.q_in_g    = 0.0
    return face
end

@testset "Isothermal enclosure with a transparent bin" begin
    T_wall = 1000.0
    κ = [0.0, 0.311, 3.462]
    cases = [(1.0, PlanckBands([5e-7, 2e-6, 6e-6, 5e-5])),
             (0.5, ConstantWeights([0.4, 0.35, 0.25]))]
    for (ε_wall, model) in cases
        mesh = RayTracingDomain2D([_transparent_face(κ, ε_wall, T_wall)], [(5, 5)])
        mesh.spectral_model = model
        mesh(1_000_000; method = :exchange, verbose = false)
        smooth!(mesh; verbose = false)

        # transparent bin: volumes absorb nothing, so their columns are exactly zero
        F1 = mesh.F_smooth[1]
        n_s = length(mesh.surface_mapping)
        @test all(F1[:, n_s+1:end] .== 0)
        @test all(isapprox.(sum(F1, dims = 2), 1.0; atol = 1e-12))

        solveEquilibrium!(mesh, mesh.F_smooth;
                          max_iters = 10_000, convergence_tol = 1e-14, verbose = false)
        err = maximum(abs(ff.T_g - T_wall) for ff in mesh.fine_mesh[1])
        @test err < 1e-9
    end
end