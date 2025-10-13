"""
Solve spectral radiative equilibrium for 3D surface-only radiation.
Key insight: View factors F are wavelength-independent, so we compute them once
and reuse for all spectral bands!
"""

function steadyStateSpectral3D!(domain::Domain3D_faces{G,P}; 
                               max_iterations::Int=500) where {G,P<:Integer}
    
    # Validate spectral setup at entry
    if isnothing(domain.wavelength_band_limits)
        error("""
        Spectral solve requires wavelength band limits to be set.
        
        Please set mesh.wavelength_band_limits before calling steadyStateSpectral3D!:
        
        Example, using logarithmic spacing:
            mesh.wavelength_band_limits = 10 .^ range(log10(0.0000001), log10(0.001), length=51)
        """)
    end
    
    # Additional validation
    if length(domain.wavelength_band_limits) < 4
        error("wavelength_band_limits must have at least 4 values (defining 3 bins)")
    end
    
    if any(domain.wavelength_band_limits .<= 0)
        error("wavelength_band_limits must all be positive (wavelengths > 0)")
    end

    if any(diff(domain.wavelength_band_limits) .<= 0)
        error("wavelength_band_limits must be strictly increasing (no duplicates)")
    end

    println("=== 3D Spectral Surface Radiation Solver ===")
    println("Spectral mode: $(domain.spectral_mode)")
    println("Number of spectral bins: $(domain.n_spectral_bins)")
    
    # Get system size
    N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    
    # KEY INSIGHT: Use the SAME F matrix for all spectral bands!
    # View factors are geometry-only and wavelength-independent
    F = domain.F
    
    println("\nComputing GERT matrices for each spectral band...")
    println("(Using same view factor matrix F for all bands)")
    
    # Build system matrices for each spectral band
    # Each band has different epsilon, so different B, K, S_infty, A, R, C, D, M
    D_matrices = Vector{Matrix{G}}()
    C_matrices = Vector{Matrix{G}}()
    M_matrices = Vector{Matrix{G}}()
    
    for bin in 1:domain.n_spectral_bins
        println("  Building matrices for spectral bin $bin...")
        B, S_infty, A, R, C, D, M = build_system_matrices3D!(domain, F; spectral_bin=bin)
        push!(D_matrices, D)
        push!(C_matrices, C)
        push!(M_matrices, M)
    end
    
    println("\nAssembling block matrix structure...")
    # Build block matrix (same structure as 2D)
    if G <: Measurement
        block_matrix = zeros(G, (domain.n_spectral_bins + 1) * N_surfs, 
                            domain.n_spectral_bins * N_surfs)
    else
        block_matrix = spzeros((domain.n_spectral_bins + 1) * N_surfs, 
                              domain.n_spectral_bins * N_surfs)
    end
    
    for i = 1:(domain.n_spectral_bins + 1)
        for j = 1:domain.n_spectral_bins
            row_start = (i - 1) * N_surfs + 1
            row_end = i * N_surfs
            col_start = (j - 1) * N_surfs + 1
            col_end = j * N_surfs
            
            if i == 1
                block_matrix[row_start:row_end, col_start:col_end] = M_matrices[j]
            elseif i == j + 1
                block_matrix[row_start:row_end, col_start:col_end] = D_matrices[j]
            end
        end
    end
    # Outside iteration loop - factorize once
    Factorization = qr(block_matrix)
    
    # Setup boundary conditions
    println("Setting up boundary conditions...")
    boundary, temperatures, emissive = setup_boundary_conditions_3D(domain)
    
    # Iteration loop
    sol_j = zeros(G, domain.n_spectral_bins * N_surfs)
    previous_sol_j = zeros(G, domain.n_spectral_bins * N_surfs)
    emit_frac = getBinsEmissionFractions_3D(domain, temperatures)
    
    record_convergence = G[]
    
    println("\nStarting spectral iteration...")
    for iter = 1:max_iterations
        
        # Update emissive powers and temperatures
        emissive = updateSpectralEmission_3D!(domain, iter, D_matrices, sol_j, 
                                             emit_frac, temperatures, emissive)
        temperatures = updateTemperaturesSpectral_3D!(domain, emissive, emit_frac)
        emit_frac = getBinsEmissionFractions_3D(domain, temperatures)
        
        # Solve block system
        emissive_pow_vec = reduce(vcat, [emissive for _ in 1:domain.n_spectral_bins])
        b_e_matrix = Diagonal(reduce(vcat, [boundary; emissive_pow_vec]))
        sol_j .= Factorization \ (b_e_matrix*[ones(G, N_surfs); emit_frac[:]])
        
        # Check convergence
        convergence_error = maximum(abs.(sol_j - previous_sol_j))
        push!(record_convergence, convergence_error)
        previous_sol_j .= sol_j
        
        println("Iteration $iter: convergence error = $convergence_error")
        
        if iter > 1 && convergence_error < 100_000*eps(Float64)
            println("Converged after $iter iterations")
            
            # Compute energy conservation error for each band
            domain.energy_error = G.([sum(C_matrices[i] * sol_j[(i - 1) * N_surfs + 1:i * N_surfs]) 
                               for i in 1:domain.n_spectral_bins])
            println("Energy conservation errors by band: $(domain.energy_error)")
            
            break
        end
        
        if iter == max_iterations
            println("Warning: Maximum iterations reached. Final error = $convergence_error")
        end
    end
    
    # AFTER convergence: Write spectral results for each band to mesh FIRST
    println("\nWriting spectral results to mesh...")
    for bin in 1:domain.n_spectral_bins
        # Extract solution for this band
        bin_start = (bin - 1) * N_surfs + 1
        bin_end = bin * N_surfs
        j_bin = sol_j[bin_start:bin_end]
        
        # Compute e, r, g_a from GERT matrices and radiosity
        e_bin = D_matrices[bin] * j_bin           # e = D * j
        r_bin = j_bin - e_bin                     # r = j - e = R' * j
        g_a_bin = j_bin - C_matrices[bin] * j_bin - r_bin  # g_a = A' * j = j - C*j - r
        
        # Write to mesh for this band
        surf_count = 0
        for superface in domain.facesMesh
            for subface in superface.subFaces
                surf_count += 1
                
                if isa(subface.e_w, Vector)
                    subface.j_w[bin] = j_bin[surf_count]
                    subface.e_w[bin] = max(e_bin[surf_count], 0.0)
                    subface.r_w[bin] = max(r_bin[surf_count], 0.0)
                    subface.g_a_w[bin] = max(g_a_bin[surf_count], 0.0)
                    subface.g_w[bin] = subface.g_a_w[bin] + subface.r_w[bin]
                    subface.i_w[bin] = j_bin[surf_count] / (π * subface.area)
                end
            end
        end
    end
    
    # THEN update scalar temperatures and heat fluxes from spectral results
    println("Computing final scalar temperatures and heat fluxes...")
    update_scalar_temperatures_3D!(domain)
    
    println("=== 3D Spectral Solution Complete ===")
end

# Helper functions (adapted from 2D versions)

function setup_boundary_conditions_3D(domain::Domain3D_faces{G,P}) where {G,P}
    
    N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    boundary = zeros(G, N_surfs)
    temperatures = zeros(G, N_surfs)
    emissive = zeros(G, N_surfs)
    
    # Get initial temperatures
    surf_count = 0
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            if subface.T_in_w > -0.1
                temperatures[surf_count] = subface.T_in_w
            end
        end
    end
    
    # Calculate emissive fractions
    emit_frac = getBinsEmissionFractions_3D(domain, temperatures)
    
    # Set boundary conditions
    surf_count = 0
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            
            epsilon_vals = isa(subface.epsilon, Vector) ? subface.epsilon : [subface.epsilon]
            weighted_epsilon = sum([epsilon_vals[i] * emit_frac[surf_count, i] 
                                   for i in 1:domain.n_spectral_bins])
            
            if subface.T_in_w > -0.1
                # Known temperature
                emissive[surf_count] = weighted_epsilon * subface.area * 
                                      STEFAN_BOLTZMANN * subface.T_in_w^4
                boundary[surf_count] = emissive[surf_count]
            else
                # Known flux
                emissive[surf_count] = weighted_epsilon * subface.area * 
                                      STEFAN_BOLTZMANN * maximum(temperatures)^4
                boundary[surf_count] = subface.q_in_w
            end
        end
    end
    
    return boundary, temperatures, emissive
end

function getBinsEmissionFractions_3D(domain::Domain3D_faces{G,P},
                                    temperatures::Vector{G}) where {G,P}
    
    N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    emit_frac = zeros(G, N_surfs, domain.n_spectral_bins)
    
    surf_count = 0
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            
            F_0_to_lambda_T_prev = 0
            for spectral_pos in 1:domain.n_spectral_bins
                F_0_to_lambda_T = emitFracBlackBody_element_spectrum(domain.wavelength_band_limits, 
                                                                     temperatures[surf_count], 
                                                                     spectral_pos)
                if spectral_pos == domain.n_spectral_bins
                    emit_frac[surf_count, spectral_pos] = 1.0 - F_0_to_lambda_T
                else
                    emit_frac[surf_count, spectral_pos] = F_0_to_lambda_T - F_0_to_lambda_T_prev
                end
                F_0_to_lambda_T_prev = F_0_to_lambda_T
            end
        end
    end
    
    return emit_frac
end

function updateSpectralEmission_3D!(domain::Domain3D_faces{G,P}, iter::Int, 
                                   D_matrices::Vector{Matrix{G}}, sol_j::Vector{G}, 
                                   emit_frac::Matrix{G}, temperatures::Vector{G}, 
                                   emissive::Vector{G}) where {G,P}
    
    N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    
    if iter > 1
        # Update from previous iteration
        Ds_combined = hcat([D_matrices[i] for i in 1:domain.n_spectral_bins]...)
        emissive .= max.(Ds_combined * sol_j, 10*eps(G))
    else
        # Initial guess - MUST use weighted epsilon!
        surf_count = 0
        for superface in domain.facesMesh
            for subface in superface.subFaces
                surf_count += 1
                
                # # Get weighted epsilon based on temperature and emission fractions
                epsilon_vals = isa(subface.epsilon, Vector) ? subface.epsilon : [subface.epsilon]
                weighted_epsilon = sum([epsilon_vals[i] * emit_frac[surf_count, i] 
                                       for i in 1:domain.n_spectral_bins])
                
                if subface.T_in_w < -0.1
                    # Known flux - use max temperature
                    emissive[surf_count] = weighted_epsilon * subface.area * STEFAN_BOLTZMANN * maximum(temperatures)^4
                else
                    # Known temperature
                    emissive[surf_count] = weighted_epsilon * subface.area * STEFAN_BOLTZMANN * subface.T_in_w^4 
                end
            end
        end
    end
    
    return emissive
end

function updateTemperaturesSpectral_3D!(domain::Domain3D_faces{G,P}, 
                                       emissive::Vector{G}, 
                                       emit_frac::Matrix{G}) where {G,P}
    
    N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    temperatures = zeros(G, N_surfs)
    
    surf_count = 0
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            
            if subface.T_in_w < -0.1
                # Radiative equilibrium - solve for scalar temperature
                epsilon_vals = isa(subface.epsilon, Vector) ? subface.epsilon : [subface.epsilon]
                weighted_epsilon = sum([epsilon_vals[j] * emit_frac[surf_count, j] 
                                       for j in 1:domain.n_spectral_bins])
                
                # Update scalar temperature
                subface.T_w = (emissive[surf_count] / 
                              (weighted_epsilon * STEFAN_BOLTZMANN * subface.area))^(1/4)
                temperatures[surf_count] = subface.T_w
            else
                # Prescribed temperature (scalar)
                subface.T_w = subface.T_in_w
                temperatures[surf_count] = subface.T_in_w
            end
        end
    end
    
    return temperatures
end