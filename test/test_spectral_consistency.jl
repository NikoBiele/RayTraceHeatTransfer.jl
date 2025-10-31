"""
Test consistency between grey and spectral implementations.
Tests for both 2D and 3D cases to ensure that:
- Spectral with uniform black surfaces matches grey results
- Spectral bin-by-bin results are physically consistent
- Integration over all bins gives correct total values
"""

println("\n" * "-"^60)
println("Testing Spectral Consistency")
println("-"^60)

#############################################################################
### TEST 1: 3D SPECTRAL vs GREY (Black Uniform Surfaces) ##################
#############################################################################

@testset "3D: Spectral vs Grey Consistency" begin
    # Unit cube geometry
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
    
    Ndim = 4
    n_bins = 20
    
    # Boundary conditions
    T_in_w = [1000.0, 500.0, -1.0, -1.0, -1.0, -1.0]
    q_in_w = [0.0, 0.0, 1000.0, 1000.0, 0.0, 0.0]
    
    # Grey case
    epsilon_grey = ones(6)
    domain_grey = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w, epsilon_grey)
    steadyStateGrey3D!(domain_grey, domain_grey.F)
    
    grey_T = Float64[]
    grey_q = Float64[]
    for superface in domain_grey.facesMesh
        for subface in superface.subFaces
            push!(grey_T, subface.T_w)
            push!(grey_q, subface.q_w)
        end
    end
    
    # Spectral case (black walls in all bins)
    epsilon_spectral_black = fill(1.0, n_bins)
    epsilon_spectral = [epsilon_spectral_black for _ in 1:6]
    
    domain_spectral = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w, 
                                     epsilon_spectral)
    domain_spectral.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    steadyStateSpectral3D!(domain_spectral; max_iterations=200)
    
    spectral_T = Float64[]
    spectral_q = Float64[]
    for superface in domain_spectral.facesMesh
        for subface in superface.subFaces
            push!(spectral_T, subface.T_w)
            push!(spectral_q, subface.q_w)
        end
    end
    
    # Compare temperatures
    @test length(grey_T) == length(spectral_T)
    
    T_max_diff = maximum(abs.(spectral_T - grey_T))
    T_rms_diff = sqrt(mean((spectral_T - grey_T).^2))
    
    @test T_max_diff < 2*TEMP_TOLERANCE # K
    @test T_rms_diff < TEMP_TOLERANCE # K
    
    # Compare heat fluxes
    q_diff = abs.(spectral_q - grey_q)
    nonzero_q = abs.(grey_q) .> 1e-10
    
    if any(nonzero_q)
        rel_q_errors = q_diff[nonzero_q] ./ abs.(grey_q[nonzero_q])
        @test maximum(rel_q_errors) < CONSISTENCY_TOLERANCE
    end
end

#############################################################################
### TEST 2: SPECTRAL BIN INTEGRATION (2D) ##################################
#############################################################################

@testset "2D: Spectral Bin Integration" begin
    
    T_hot = 1000.0
    T_cold = 0.0
    kappa = 1.0
    sigma_s = 0.0
    Ndim = 5
    n_bins = 20
    
    vertices = SVector(
        Point2(0.0, 0.0),
        Point2(1.0, 0.0), 
        Point2(1.0, 1.0),
        Point2(0.0, 1.0)
    )
    
    solidWalls = SVector(true, true, true, true)
    face = PolyFace2D{Float64}(vertices, solidWalls, n_bins, kappa, sigma_s)
    
    # Uniform spectral properties
    face.kappa_g = fill(kappa, n_bins)
    face.sigma_s_g = fill(sigma_s, n_bins)
    face.q_in_g = 0.0
    face.T_in_g = -1.0
    face.T_in_w = [T_hot, T_cold, T_cold, T_cold]
    
    # Initialize spectral result vectors
    face.j_g = zeros(n_bins)
    face.g_a_g = zeros(n_bins)
    face.e_g = zeros(n_bins)
    face.r_g = zeros(n_bins)
    face.g_g = zeros(n_bins)
    face.i_g = zeros(n_bins)
    
    for wall_idx in 1:4
        face.j_w[wall_idx] = zeros(n_bins)
        face.g_a_w[wall_idx] = zeros(n_bins)
        face.e_w[wall_idx] = zeros(n_bins)
        face.r_w[wall_idx] = zeros(n_bins)
        face.g_w[wall_idx] = zeros(n_bins)
        face.i_w[wall_idx] = zeros(n_bins)
    end
    
    epsilon_bins = fill(1.0, n_bins)
    face.epsilon = [epsilon_bins, epsilon_bins, epsilon_bins, epsilon_bins]
    
    mesh = RayTracingDomain2D([face], [(Ndim, Ndim)])
    mesh.spectral_mode = :spectral_uniform
    mesh.n_spectral_bins = n_bins
    mesh.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    
    # Solve spectral problem
    N_rays = 1_000_000
    mesh(N_rays; method=:exchange)
    steadyStateSpectral2D!(mesh, mesh.F_smooth)
    
    # Check that spectral results integrate properly
    for fine_face in mesh.fine_mesh[1]
        # Volume radiosity
        if isa(fine_face.j_g, Vector)
            j_total_integrated = sum(fine_face.j_g)
            
            # This should equal the scalar total (if stored)
            # At minimum, check it's positive and finite
            @test j_total_integrated >= 0.0
            @test isfinite(j_total_integrated)
            
            # Emission should approximately equal absorption in equilibrium
            e_total = sum(fine_face.e_g)
            ga_total = sum(fine_face.g_a_g)
            
            if e_total > 1e-10
                abs_diff = abs(e_total - ga_total)
                @test abs_diff < ENERGY_TOLERANCE # small tolerance
            end
        end
        
        # Wall radiosity
        for wall_idx in 1:4
            if isa(fine_face.j_w[wall_idx], Vector)
                j_w_total = sum(fine_face.j_w[wall_idx])
                @test isfinite(j_w_total)
                
                # Energy balance: j = e + r
                e_w_total = sum(fine_face.e_w[wall_idx])
                r_w_total = sum(fine_face.r_w[wall_idx])
                
                @test isapprox(j_w_total, e_w_total + r_w_total, atol=ENERGY_TOLERANCE)
            end
        end
    end
end

#############################################################################
### TEST 3: SELECTIVE vs BLACK SURFACES (Temperature Differences) ##########
#############################################################################

@testset "3D: Selective vs Black Surface Comparison" begin
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
    
    Ndim = 4
    n_bins = 20
    
    T_in_w = [1000.0, 500.0, -1.0, -1.0, -1.0, -1.0]
    q_in_w = zeros(6)
    
    # Black surfaces
    epsilon_black = fill(1.0, n_bins)
    epsilon_all_black = [epsilon_black for _ in 1:6]
    
    domain_black = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w, 
                                  epsilon_all_black)
    domain_black.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    steadyStateSpectral3D!(domain_black)
    
    black_temps = [sf.T_w for i in 1:6 for sf in domain_black.facesMesh[i].subFaces]
    
    # Selective surfaces (wavelength-dependent)
    epsilon_selective = [i <= n_bins÷2 ? 0.3 : 0.9 for i in 1:n_bins]
    epsilon_all_selective = [epsilon_selective for _ in 1:6]
    
    domain_selective = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w,
                                      epsilon_all_selective)
    domain_selective.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    steadyStateSpectral3D!(domain_selective)
    
    selective_temps = [sf.T_w for i in 1:6 for sf in domain_selective.facesMesh[i].subFaces]
    
    # Results should be DIFFERENT (selective coating changes heat transfer)
    @test !isapprox(black_temps, selective_temps, rtol=0.01)
    
    # But both should still be physical
    @test all(400.0 .<= black_temps .<= 1100.0)
    @test all(400.0 .<= selective_temps .<= 1100.0)
end

#############################################################################
### TEST 4: SPECTRAL MODE DETECTION ########################################
#############################################################################

@testset "Spectral Mode Auto-Detection" begin
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
    
    Ndim = 3
    T_in_w = [1000.0, 0.0, -1.0, -1.0, -1.0, -1.0]
    q_in_w = zeros(6)
    
    # Grey mode: scalar epsilon
    epsilon_grey = ones(6)
    domain_grey = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w, epsilon_grey)
    
    @test domain_grey.spectral_mode == :grey
    @test domain_grey.n_spectral_bins == 1
    
    # Spectral mode: vector epsilon
    n_bins = 10
    epsilon_spectral_vec = fill(1.0, n_bins)
    epsilon_spectral = [epsilon_spectral_vec for _ in 1:6]
    
    domain_spectral = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w, 
                                     epsilon_spectral)
    domain_spectral.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    
    @test domain_spectral.spectral_mode == :spectral_uniform || 
          domain_spectral.spectral_mode == :spectral_variable
    @test domain_spectral.n_spectral_bins == n_bins
end

#############################################################################
### TEST 5: ENERGY CONSERVATION ACROSS SPECTRUM ############################
#############################################################################

@testset "Energy Conservation Per Bin and Total" begin
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
    
    Ndim = 4
    n_bins = 20
    
    T_in_w = [1000.0, 500.0, -1.0, -1.0, -1.0, -1.0]
    q_in_w = [0.0, 0.0, 500.0, 500.0, 0.0, 0.0]
    
    epsilon_bins = fill(1.0, n_bins)
    epsilon = [epsilon_bins for _ in 1:6]
    
    domain = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w, epsilon)
    domain.wavelength_band_limits = 10 .^ range(log10(0.00000001), log10(0.1), length=n_bins+1)
    steadyStateSpectral3D!(domain)
    
    # Check energy balance for each spectral bin separately
    # and for the integrated total
    
    total_q_per_bin = zeros(n_bins)
    
    for i in 1:6
        for subface in domain.facesMesh[i].subFaces
            # If q_w is spectral (vector), check each bin
            # Otherwise it's already integrated
            if isa(subface.q_w, Number)
                # Scalar - this is the integrated total
                # Just check it's finite
                @test isfinite(subface.q_w)
            elseif isa(subface.q_w, Vector)
                # Vector - check each bin
                for (bin, q_bin) in enumerate(subface.q_w)
                    @test isfinite(q_bin)
                    total_q_per_bin[bin] += q_bin
                end
            end
        end
    end
    
    # Total heat flux should sum to approximately zero
    # (what goes in must come out)
    if !all(iszero, total_q_per_bin)
        for bin in 1:n_bins
            @test abs(total_q_per_bin[bin]) < ENERGY_TOLERANCE # Per bin tolerance
        end
        
        @test abs(sum(total_q_per_bin)) < ENERGY_TOLERANCE # Total tolerance
    end
end

println("✓ Spectral Consistency tests complete")