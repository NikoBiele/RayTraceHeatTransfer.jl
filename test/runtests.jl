using RayTraceHeatTransfer
using Test

@testset "RayTraceHeatTransfer.jl" begin
    println("\n" * "="^80)
    println("STARTING TEST SUITE")
    println("="^80)

    # Test 3D transparent surfaces (view factors and heat transfer)
    @testset "3D Surfaces (Transparent)" begin
        include("test_3d_viewfactors.jl")
        include("test_3d_heat_transfer.jl")
    end
    
    # Test 2D participating media
    @testset "2D Participating Media" begin
        include("test_2d_grey.jl")
        include("test_2d_grey_reflecting.jl")
        include("test_2d_spectral.jl")
    end
    
    # Test spectral capabilities
    @testset "Spectral Radiation" begin
        include("test_spectral_consistency.jl")
    end

    # Test 2d diffusion (sparse/dense x grey/spectral)
    @testset "2D Diffusion" begin
        include("test_2d_diffusion.jl")
    end

    # Test 2d spectral (sparse/dense)
    @testset "2D spectral sparse/dense" begin
        include("test_2d_spectral_dense_sparse.jl")
    end

    @testset "2D triangular mesh" begin
        include("test_triangle_mesh.jl")
    end

    @testset "Shared F grey vs spectral" begin
        include("test_shared_F_grey_vs_spectral.jl")
    end

    @testset "Selective emissivity decoupling" begin
        include("test_selective_emissivity_decoupling.jl")
    end

    @testset "3D mixed epsilon grey vs spectral" begin
        include("test_3d_mixed_eps_grey_vs_spectral.jl")
    end

    @testset "Spectral isothermal cavity" begin
        include("test_spectral_isothermal.jl")
    end

    @testset "Spectral bin refinement" begin
        include("test_spectral_bin_refinement.jl")
    end

    @testset "Geometry" begin
        include("test_geometry.jl")
    end

    @testset "Reproducible" begin
        include("test_reproducible.jl")
    end

    @testset "Transparent bin" begin
        include("test_transparent_bin.jl")
    end

    @testset "PiecewiseBands" begin
        include("test_piecewise_bands.jl")
    end

    @testset "adaptiveSpectralBins" begin
        include("test_adaptive_spectral_bins.jl")
    end

    @testset "Pathlength tracer" begin
        include("test_pathlength_tracer.jl")
    end
    
    @testset "Line-by-line solve validation" begin
        include("test_lbl_slab.jl")
    end

    @testset "Directional models" begin
        include("test_directional_models.jl")
    end

    @testset "Directional solvers" begin
        include("test_directional_solvers.jl")
    end

    @testset "Global solution field" begin
        include("test_solution_field.jl")
    end

    @testset "Anisotropic slab validation" begin
        include("test_anisotropic_slab.jl")
    end

    println("\n" * "="^80)
    println("TEST SUITE COMPLETE")
    println("="^80)
end