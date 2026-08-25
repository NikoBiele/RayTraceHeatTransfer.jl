println("\n" * "-"^60)
println("Testing 2D triangular meshing (circular meshes)")
println("-"^60)

using RayTraceHeatTransfer
using Test
using LinearAlgebra
using StatsBase
using StaticArrays
using GeometryBasics
using SparseArrays
using ConvolutionInterpolations

# Circle from triangles: N_seg wedges sharing a center vertex.
# Wedge j: (center, P_j, P_{j+1}); rim edge solid, both spokes open.

N_seg  = 16;
R      = 1.0;
T_hot  = 1000.0;
kappa  = 1.0;

rim(j) = Point2(R * cos(2π * (j - 1) / N_seg), R * sin(2π * (j - 1) / N_seg));

@testset "Isothermal circular triangle mesh" begin
    # ---- Case 1: isothermal rim -> gas equilibrium T must equal T_hot ----------
    faces     = PolyVolume2D{Float64}[];
    divisions = Tuple{Int,Int}[];
    for j in 1:N_seg
        verts = SVector(Point2(0.0, 0.0), rim(j), rim(j + 1))
        solidwalls = SVector(false, true, false)      # spoke, rim, spoke
        face = PolyVolume2D{Float64}(verts, solidwalls, 1, kappa, 0.0)
        face.T_in_w  = [0.0, T_hot, 0.0]              # only rim wall is real
        face.epsilon = [1.0, 1.0, 1.0]
        face.T_in_g  = -1.0
        face.q_in_g  = 0.0
        push!(faces, face)
        push!(divisions, (2, 2))
    end;
    mesh = RayTracingDomain2D(faces, divisions);
    mesh(10_000_000; method = :exchange);
    smooth!(mesh)
    solveEquilibrium!(mesh, mesh.F_smooth);
    T_g = [ff.T_g for fine in mesh.fine_mesh for ff in fine];
    # test for isothermal
    @test extrema(T_g)[1] > T_hot - 1e-3
    @test extrema(T_g)[2] < T_hot + 1e-3
end

@testset "Center temperature and symmetry of circular triangular mesh" begin
   # ---- Case 2: half rim hot, half cold -> center symmetry limit --------------
    faces2     = PolyVolume2D{Float64}[];
    divisions2 = Tuple{Int,Int}[];
    for j in 1:N_seg
        verts = SVector(Point2(0.0, 0.0), rim(j), rim(j + 1))
        face = PolyVolume2D{Float64}(verts, SVector(false, true, false), 1, kappa, 0.0)
        face.T_in_w  = [0.0, j <= N_seg ÷ 2 ? T_hot : 0.0, 0.0]
        face.epsilon = [1.0, 1.0, 1.0]
        face.T_in_g  = -1.0
        face.q_in_g  = 0.0
        push!(faces2, face)
        push!(divisions2, (11, 11))
    end;
    mesh2 = RayTracingDomain2D(faces2, divisions2);
    mesh2(10_000_000; method = :exchange);
    smooth!(mesh2)
    solveEquilibrium!(mesh2, mesh2.F_smooth);

    T_limit = ((T_hot^4 + 0.0^4) / 2)^(1/4);      # ≈ 840.896 K
    # gas elements adjacent to the center
    T_g_mid = [fine[1].T_g for fine in mesh2.fine_mesh];
    @test abs(T_limit - mean(T_g_mid)) < 2.0

    # test symmetry
    T_g_all = [ff.T_g for fine in mesh2.fine_mesh for ff in fine];
    @test abs(mean((T_g_all ./ T_hot).^4) - 0.5) < 0.1
     
end