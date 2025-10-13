"""
Test 2D spectral participating media ray tracing.
Tests include:
- Spectral uniform (constant properties across wavelengths)
- Spectral variable (wavelength-dependent properties)
- Comparison with grey case when appropriate
- Method consistency (exchange factor vs direct)
"""

println("\n" * "-"^60)
println("Testing 2D Spectral Participating Media")
println("-"^60)

#############################################################################
### HELPER FUNCTIONS #######################################################
#############################################################################

function create_spectral_uniform_mesh(; T_hot=1000.0, T_cold=0.0, kappa=1.0, 
                                       sigma_s=0.0, epsilon_bins=nothing, 
                                       Ndim=5, n_bins=10)
    """
    Create a 2D square with uniform spectral properties.
    """
    
    vertices = SVector(
        Point2(0.0, 0.0),
        Point2(1.0, 0.0), 
        Point2(1.0, 1.0),
        Point2(0.0, 1.0)
    )
    
    solidWalls = SVector(true, true, true, true)
    face = PolyFace2D{Float64}(vertices, solidWalls, n_bins, kappa, sigma_s)
    
    # Set spectral properties (uniform across bins)
    face.kappa_g = fill(kappa, n_bins)
    face.sigma_s_g = fill(sigma_s, n_bins)
    face.q_in_g = 0.0
    face.T_in_g = -1.0
    face.T_in_w = [T_hot, T_cold, T_cold, T_cold]
    
    # Convert result fields to vectors
    face.j_g = zeros(n_bins)
    face.g_a_g = zeros(n_bins)
    face.e_g = zeros(n_bins)
    face.r_g = zeros(n_bins)
    face.g_g = zeros(n_bins)
    face.i_g = zeros(n_bins)
    
    # Wall fields to vectors
    for wall_idx in 1:4
        face.j_w[wall_idx] = zeros(n_bins)
        face.g_a_w[wall_idx] = zeros(n_bins)
        face.e_w[wall_idx] = zeros(n_bins)
        face.r_w[wall_idx] = zeros(n_bins)
        face.g_w[wall_idx] = zeros(n_bins)
        face.i_w[wall_idx] = zeros(n_bins)
    end
    
    if epsilon_bins === nothing
        epsilon_bins = fill(1.0, n_bins)
    end
    face.epsilon = [epsilon_bins, epsilon_bins, epsilon_bins, epsilon_bins]
    
    mesh = RayTracingMeshOptim([face], [(Ndim, Ndim)])
    mesh.spectral_mode = :spectral_uniform
    mesh.n_spectral_bins = n_bins
    
    return mesh
end

function create_spectral_variable_mesh(; T_hot=1000.0, T_cold=0.0, base_kappa=1.0,
                                        base_sigma_s=0.0, epsilon_bins=nothing,
                                        Ndim=5, n_bins=10)
    """
    Create a 2D square with variable spectral properties.
    Extinction varies slightly across bins.
    """
    
    vertices = SVector(
        Point2(0.0, 0.0),
        Point2(1.0, 0.0), 
        Point2(1.0, 1.0),
        Point2(0.0, 1.0)
    )
    
    solidWalls = SVector(true, true, true, true)
    face = PolyFace2D{Float64}(vertices, solidWalls, n_bins, base_kappa, base_sigma_s)
    
    # Variable extinction (small variation for testing)
    face.kappa_g = [base_kappa * (1.0 + 0.01 * (i-1)/(n_bins-1)) for i in 1:n_bins]
    face.sigma_s_g = [base_sigma_s * (1.0 + 0.01 * (i-1)/(n_bins-1)) for i in 1:n_bins]
    
    face.q_in_g = 0.0
    face.T_in_g = -1.0
    face.T_in_w = [T_hot, T_cold, T_cold, T_cold]
    
    # Convert result fields to vectors
    face.j_g = zeros(n_bins)
    face.g_a_g = zeros(n_bins)
    face.e_g = zeros(n_bins)
    face.r_g = zeros(n_bins)
    face.g_g = zeros(n_bins)
    face.i_g = zeros(n_bins)
    
    # Wall fields to vectors
    for wall_idx in 1:4
        face.j_w[wall_idx] = zeros(n_bins)
        face.g_a_w[wall_idx] = zeros(n_bins)
        face.e_w[wall_idx] = zeros(n_bins)
        face.r_w[wall_idx] = zeros(n_bins)
        face.g_w[wall_idx] = zeros(n_bins)
        face.i_w[wall_idx] = zeros(n_bins)
    end
    
    if epsilon_bins === nothing
        epsilon_bins = fill(1.0, n_bins)
    end
    face.epsilon = [epsilon_bins, epsilon_bins, epsilon_bins, epsilon_bins]
    
    mesh = RayTracingMeshOptim([face], [(Ndim, Ndim)])
    mesh.spectral_mode = :spectral_variable
    mesh.n_spectral_bins = n_bins
    
    return mesh
end

#############################################################################
### TEST 1: SPECTRAL UNIFORM WITH BLACK WALLS (should match grey) ##########
#############################################################################

@testset "Spectral Uniform vs Grey (Black Walls)" begin
    T_hot = 1000.0
    T_cold = 0.0
    kappa = 1.0
    sigma_s = 0.0
    Ndim = 5
    n_bins = 10
    N_rays = 1_000_000
    
    # Grey case
    vertices = SVector(
        Point2(0.0, 0.0),
        Point2(1.0, 0.0), 
        Point2(1.0, 1.0),
        Point2(0.0, 1.0)
    )
    solidWalls = SVector(true, true, true, true)
    face_grey = PolyFace2D{Float64}(vertices, solidWalls, 1, kappa, sigma_s)
    face_grey.T_in_w = [T_hot, T_cold, T_cold, T_cold]
    face_grey.epsilon = [1.0, 1.0, 1.0, 1.0]
    face_grey.T_in_g = -1.0
    
    mesh_grey = RayTracingMeshOptim([face_grey], [(Ndim, Ndim)])
    mesh_grey(N_rays; method=:exchange)
    steadyStateGrey2D!(mesh_grey, mesh_grey.F_smooth)
    
    grey_temps = [fine_face.T_g for fine_face in mesh_grey.fine_mesh[1]]
    
    # Spectral uniform case (black walls in all bins)
    epsilon_bins = fill(1.0, n_bins)
    mesh_spectral = create_spectral_uniform_mesh(T_hot=T_hot, T_cold=T_cold,
                                                 kappa=kappa, sigma_s=sigma_s,
                                                 epsilon_bins=epsilon_bins,
                                                 Ndim=Ndim, n_bins=n_bins)
    mesh_spectral.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    mesh_spectral(N_rays; method=:exchange)
    steadyStateSpectral2D!(mesh_spectral, mesh_spectral.F_smooth)
    
    spectral_temps = [fine_face.T_g for fine_face in mesh_spectral.fine_mesh[1]]
    
    # Temperatures should match closely
    @test length(grey_temps) == length(spectral_temps)
    
    for (i, (T_grey, T_spec)) in enumerate(zip(grey_temps, spectral_temps))
        rel_diff = abs(T_grey - T_spec) / max(T_grey, 1.0)
        @test rel_diff < SPECTRAL_TOLERANCE
    end
    
    # Mean temperature comparison
    @test isapprox(mean(grey_temps), mean(spectral_temps), rtol=SPECTRAL_TOLERANCE)
end

#############################################################################
### TEST 2: EXCHANGE FACTOR VS DIRECT METHOD (Spectral Uniform) ############
#############################################################################

@testset "Method Consistency - Spectral Uniform" begin
    T_hot = 1000.0
    T_cold = 0.0
    kappa = 1.0
    sigma_s = 0.0
    Ndim = 5
    n_bins = 10
    N_rays = 1_000_000
    
    epsilon_bins = fill(1.0, n_bins)
    
    # Exchange factor method
    mesh_exchange = create_spectral_uniform_mesh(T_hot=T_hot, T_cold=T_cold,
                                                 kappa=kappa, sigma_s=sigma_s,
                                                 epsilon_bins=epsilon_bins,
                                                 Ndim=Ndim, n_bins=n_bins)
    mesh_exchange.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    mesh_exchange(N_rays; method=:exchange)
    steadyStateSpectral2D!(mesh_exchange, mesh_exchange.F_smooth)
    
    exchange_temps = [fine_face.T_g for fine_face in mesh_exchange.fine_mesh[1]]
    
    # Direct method
    mesh_direct = create_spectral_uniform_mesh(T_hot=T_hot, T_cold=T_cold,
                                               kappa=kappa, sigma_s=sigma_s,
                                               epsilon_bins=epsilon_bins,
                                               Ndim=Ndim, n_bins=n_bins)
    mesh_direct.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    mesh_direct(N_rays; method=:direct)
    
    direct_temps = [fine_face.T_g for fine_face in mesh_direct.fine_mesh[1]]
    
    # Methods should agree reasonably well
    for (T_ex, T_dir) in zip(exchange_temps, direct_temps)
        rel_diff = abs(T_ex - T_dir) / max(T_ex, 1.0)
        @test rel_diff < SPECTRAL_TOLERANCE
    end
end

#############################################################################
### TEST 3: EXCHANGE FACTOR VS DIRECT METHOD (Spectral Variable) ###########
#############################################################################

@testset "Method Consistency - Spectral Variable" begin
    T_hot = 1000.0
    T_cold = 0.0
    base_kappa = 1.0
    base_sigma_s = 0.0
    Ndim = 5
    n_bins = 10
    N_rays = 1_000_000
    
    epsilon_bins = fill(1.0, n_bins)
    
    # Exchange factor method
    mesh_exchange = create_spectral_variable_mesh(T_hot=T_hot, T_cold=T_cold,
                                                  base_kappa=base_kappa,
                                                  base_sigma_s=base_sigma_s,
                                                  epsilon_bins=epsilon_bins,
                                                  Ndim=Ndim, n_bins=n_bins)
    mesh_exchange.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    mesh_exchange(N_rays; method=:exchange)
    steadyStateSpectral2D!(mesh_exchange, mesh_exchange.F_smooth)
    
    exchange_temps = [fine_face.T_g for fine_face in mesh_exchange.fine_mesh[1]]
    
    # Direct method
    mesh_direct = create_spectral_variable_mesh(T_hot=T_hot, T_cold=T_cold,
                                                base_kappa=base_kappa,
                                                base_sigma_s=base_sigma_s,
                                                epsilon_bins=epsilon_bins,
                                                Ndim=Ndim, n_bins=n_bins)
    mesh_direct.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    mesh_direct(N_rays; method=:direct)
    
    direct_temps = [fine_face.T_g for fine_face in mesh_direct.fine_mesh[1]]
    
    # Methods should agree reasonably well
    for (T_ex, T_dir) in zip(exchange_temps, direct_temps)
        rel_diff = abs(T_ex - T_dir) / max(T_ex, 1.0)
        @test rel_diff < SPECTRAL_TOLERANCE
    end
end

#############################################################################
### TEST 4: SPECTRAL ENERGY CONSERVATION ###################################
#############################################################################

@testset "Spectral Energy Conservation" begin
    T_hot = 1000.0
    T_cold = 0.0
    kappa = 1.0
    sigma_s = 0.0
    Ndim = 5
    n_bins = 10
    N_rays = 1_000_000
    
    epsilon_bins = fill(1.0, n_bins)
    
    mesh = create_spectral_uniform_mesh(T_hot=T_hot, T_cold=T_cold,
                                       kappa=kappa, sigma_s=sigma_s,
                                       epsilon_bins=epsilon_bins,
                                       Ndim=Ndim, n_bins=n_bins)
    mesh.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    mesh(N_rays; method=:exchange)
    steadyStateSpectral2D!(mesh, mesh.F_smooth)
    
    # Check energy error for each bin if available
    if !isnothing(mesh.energy_error)
        if isa(mesh.energy_error, Vector)
            for (bin, err) in enumerate(mesh.energy_error)
                @test abs(err) < ENERGY_TOLERANCE  # Small absolute error per bin
            end
            
            # Total error across all bins
            @test abs(sum(mesh.energy_error)) < ENERGY_TOLERANCE
        else
            @test abs(mesh.energy_error) < ENERGY_TOLERANCE
        end
    end
end

#############################################################################
### TEST 5: SELECTIVE SURFACE (wavelength-dependent epsilon) ###############
#############################################################################

@testset "Selective Surface Properties" begin
    T_hot = 1000.0
    T_cold = 0.0
    kappa = 1.0
    sigma_s = 0.0
    Ndim = 4
    n_bins = 10
    N_rays = 1_000_000
    
    # Create selective emissivity (low in some bins, high in others)
    epsilon_selective = [i <= n_bins÷2 ? 0.3 : 0.9 for i in 1:n_bins]
    
    mesh = create_spectral_uniform_mesh(T_hot=T_hot, T_cold=T_cold,
                                       kappa=kappa, sigma_s=sigma_s,
                                       epsilon_bins=epsilon_selective,
                                       Ndim=Ndim, n_bins=n_bins)
    mesh.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    mesh(N_rays; method=:exchange)
    steadyStateSpectral2D!(mesh, mesh.F_smooth)
    
    # Solution should exist and be physical
    for fine_face in mesh.fine_mesh[1]
        @test 0.0 <= fine_face.T_g <= T_hot * (1.0 + ANALYTICAL_TOLERANCE)
        @test isfinite(fine_face.T_g)
        
        # Check that spectral distribution makes sense
        if isa(fine_face.j_w, Vector) && length(fine_face.j_w) > 0
            if isa(fine_face.j_w[1], Vector)
                for wall_idx in 1:4
                    # Total radiosity should be positive
                    @test sum(fine_face.j_w[wall_idx]) >= 0.0
                end
            end
        end
    end
end

#############################################################################
### TEST 6: NUMBER OF BINS VARIATION #######################################
#############################################################################

@testset "Varying Number of Spectral Bins" begin
    T_hot = 1000.0
    T_cold = 0.0
    kappa = 1.0
    sigma_s = 0.0
    Ndim = 4
    N_rays = 1_000_000
    
    bin_counts = [5, 10, 20, 50]
    mean_temps = Float64[]
    
    for n_bins in bin_counts
        @testset "n_bins = $n_bins" begin
            epsilon_bins = fill(1.0, n_bins)
            
            mesh = create_spectral_uniform_mesh(T_hot=T_hot, T_cold=T_cold,
                                               kappa=kappa, sigma_s=sigma_s,
                                               epsilon_bins=epsilon_bins,
                                               Ndim=Ndim, n_bins=n_bins)
            mesh.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
            mesh(N_rays; method=:exchange)
            steadyStateSpectral2D!(mesh, mesh.F_smooth)
            
            temps = [fine_face.T_g for fine_face in mesh.fine_mesh[1]]

            # Physical checks
            @test all(0.0 .<= trunc.(Int, temps) .<= T_hot * 1.1)
            @test all(isfinite.(temps))
        end
    end
    
end

println("✓ 2D Spectral Participating Media tests complete")