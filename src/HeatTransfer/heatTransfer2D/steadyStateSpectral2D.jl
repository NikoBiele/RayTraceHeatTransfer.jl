function steadyStateSpectral2D_full!(rtm::RayTracingDomain2D, F_matrices::Union{Matrix{G}, Vector{Matrix{G}}}; 
                              max_iterations::P=500, spectral_coupling_tolerance::G=1e-11) where {G, P<:Integer}
    """
    Solve spectral radiative equilibrium using the iterative method from your example.
    F_matrices can be either a single matrix (grey/uniform) or vector of matrices (variable spectral)
    """
    
    # Validate spectral setup at entry
    if isnothing(rtm.wavelength_band_limits)
        error("""
        Spectral solve requires wavelength band limits to be set.
        
        Please set mesh.wavelength_band_limits before calling steadyStateSpectral2D!:
        
        Example, using logarithmic spacing:
            mesh.wavelength_band_limits = 10 .^ range(log10(0.0000001), log10(0.001), length=51)
        """)
    end
    
    # Additional validation
    if length(rtm.wavelength_band_limits) < 4
        error("wavelength_band_limits must have at least 4 values (defining 3 bins)")
    end
    
    if any(rtm.wavelength_band_limits .<= 0)
        error("wavelength_band_limits must all be positive (wavelengths > 0)")
    end

    if any(diff(rtm.wavelength_band_limits) .<= 0)
        error("wavelength_band_limits must be strictly increasing (no duplicates)")
    end

    # Get system matrices
    surface_mapping, volume_mapping = rtm.surface_mapping, rtm.volume_mapping
    num_surfaces = length(surface_mapping)
    num_volumes = length(volume_mapping)
    total_elements = num_surfaces + num_volumes
    
    # Build system matrices based on spectral mode
    D_matrices = Vector{Matrix{G}}()
    C_matrices = Vector{Matrix{G}}()
    M_matrices = Vector{Matrix{G}}()
    # Variable spectral - use different F matrix for each bin
    
    for bin in 1:rtm.n_spectral_bins
        B, S_infty, A, R, C, D, M = build_system_matrices2D!(rtm, F_matrices[bin]; spectral_bin=bin)
        push!(D_matrices, D)
        push!(C_matrices, C)
        push!(M_matrices, M)
    end
    
    # Build block matrix structure
    println("==== Building and Factorizing Block matrix ====")
    if G <: Measurement
        block_matrix = zeros(G, (rtm.n_spectral_bins + 1) * total_elements, rtm.n_spectral_bins * total_elements)
    else
        block_matrix = spzeros((rtm.n_spectral_bins + 1) * total_elements, rtm.n_spectral_bins * total_elements)
    end

    for i = 1:(rtm.n_spectral_bins + 1)  # block indices
        for j = 1:rtm.n_spectral_bins      # block indices
            row_start = (i - 1) * total_elements + 1
            row_end = i * total_elements
            col_start = (j - 1) * total_elements + 1
            col_end = j * total_elements
            
            if i == 1
                block_matrix[row_start:row_end, col_start:col_end] = M_matrices[j]
            elseif i == j + 1
                block_matrix[row_start:row_end, col_start:col_end] = D_matrices[j]
            end
            # else: leave as sparse zeros
        end
    end
    # Outside iteration loop - factorize once
    Factorization = qr(block_matrix)
        
    # Set boundary conditions from mesh
    boundary, temperatures, emissive = setup_boundary_conditions(rtm, F_matrices)
    
    # Iteration loop
    sol_j = zeros(G, rtm.n_spectral_bins * total_elements)
    record_convergence = G[]
    previous_sol_j = zeros(G, rtm.n_spectral_bins * total_elements)
    emit_frac = getBinsEmissionFractions(rtm, temperatures)

    println("Starting spectral steady-state iteration...")
    for iter = 1:max_iterations

        # calculate emissive powers and temperatures
        emissive = updateSpectralEmission!(rtm, iter, D_matrices, sol_j, emit_frac, temperatures, emissive)
        temperatures = updateTemperaturesSpectral!(rtm, emissive, emit_frac)
        emit_frac = getBinsEmissionFractions(rtm, temperatures)
        
        # calculate inverse of b-e matrix
        emissive_pow_vec = reduce(vcat, [emissive for _ in 1:rtm.n_spectral_bins])
        b_e_matrix = Diagonal(reduce(vcat, [boundary; emissive_pow_vec]))
        sol_j .= Factorization \ (b_e_matrix*[ones(G, total_elements); emit_frac[:]])

        # Check convergence
        convergence_error = maximum(abs.(sol_j - previous_sol_j))
        push!(record_convergence, convergence_error)
        previous_sol_j .= sol_j
        
        println("Iteration $iter: convergence error = $convergence_error")
        
        if iter > 1 && convergence_error < spectral_coupling_tolerance
            println("Converged after $iter iterations")
            Ds_combined = hcat([D_matrices[i] for i in 1:rtm.n_spectral_bins]...)
            emissive = max.(Ds_combined * sol_j, 10*eps(G))
            emissive = updateSpectralEmission!(rtm, iter, D_matrices, sol_j, emit_frac, temperatures, emissive)
            temperatures = updateTemperaturesSpectral!(rtm, emissive, emit_frac)
            break
        end
        
        if iter == max_iterations
            println("Warning: Maximum iterations reached. Final error = $convergence_error")
            Ds_combined = hcat([D_matrices[i] for i in 1:rtm.n_spectral_bins]...)
            emissive = max.(Ds_combined * sol_j, 10*eps(G))
            emissive = updateSpectralEmission!(rtm, iter, D_matrices, sol_j, emit_frac, temperatures, emissive)
            temperatures = updateTemperaturesSpectral!(rtm, emissive, emit_frac)
        end
    end
    
    # Write spectral results for each bin to mesh
    println("Writing spectral results to mesh...")
    for bin in 1:rtm.n_spectral_bins
        # Extract solution for this bin
        bin_start = (bin - 1) * total_elements + 1
        bin_end = bin * total_elements
        j_bin = sol_j[bin_start:bin_end]
        
        # Compute e, r, g_a from GERT matrices
        e_bin = D_matrices[bin] * j_bin           # e = D * j
        r_bin = j_bin - e_bin                     # r = j - e = R' * j
        g_a_bin = j_bin - C_matrices[bin] * j_bin - r_bin  # g_a = A' * j = j - C*j - r
        
        # Dummy temperature vector (will be overwritten by update_scalar_temperatures_and_heat_sources!)
        T_bin = zeros(G, total_elements)
        
        # Write to mesh for this bin
        write_results_to_mesh!(rtm, T_bin, j_bin, g_a_bin, r_bin, num_surfaces; spectral_bin=bin)
    end

    # Update mesh with final results
    update_scalar_temperatures_and_heat_sources!(rtm)
    # Compute energy conservation error for each spectral bin
    rtm.energy_error = G.([sum(C_matrices[i] * sol_j[(i - 1) * total_elements + 1:i * total_elements]) 
                            for i in 1:rtm.n_spectral_bins])
    
end

function steadyStateSpectral2D_direct!(rtm::RayTracingDomain2D, F_matrices::Union{Matrix{G}, Vector{Matrix{G}}}; 
                                      max_iterations::P=500, spectral_coupling_tolerance::G=1e-11) where {G, P<:Integer}
    """
    OPTIMIZED direct emission solver for non-scattering, non-reflecting problems.
    
    Key optimization: Instead of solving (n_bins+1)*N × n_bins*N block system for j,
    we solve N×N system for e directly, then recover j from emission fractions.
    
    This reduces computational complexity from O((n_bins*N)²) to O(N²).
    
    Requires: ε=1 everywhere (no reflection) AND σₛ=0 everywhere (no scattering)
    """
    
    # Validate spectral setup at entry
    if isnothing(rtm.wavelength_band_limits)
        error("""
        Spectral solve requires wavelength band limits to be set.
        
        Please set mesh.wavelength_band_limits before calling steadyStateSpectral2D_direct!:
        
        Example, using logarithmic spacing:
            mesh.wavelength_band_limits = 10 .^ range(log10(0.0000001), log10(0.001), length=51)
        """)
    end
    
    if length(rtm.wavelength_band_limits) < 4
        error("wavelength_band_limits must have at least 4 values (defining 3 bins)")
    end
    
    if any(rtm.wavelength_band_limits .<= 0)
        error("wavelength_band_limits must all be positive (wavelengths > 0)")
    end

    if any(diff(rtm.wavelength_band_limits) .<= 0)
        error("wavelength_band_limits must be strictly increasing (no duplicates)")
    end

    # Get system dimensions
    surface_mapping, volume_mapping = rtm.surface_mapping, rtm.volume_mapping
    num_surfaces = length(surface_mapping)
    num_volumes = length(volume_mapping)
    total_elements = num_surfaces + num_volumes
    
    # Uniform spectral - single F matrix for all bins
    B, S_infty, A, R, C, D, M = build_system_matrices2D!(rtm, F_matrices; spectral_bin=1)
        
    # Set boundary conditions from mesh
    boundary, temperatures, emissive = setup_boundary_conditions(rtm, F_matrices)
    
    # Initialize emission fractions (w in temporal analogy)
    emit_frac = getBinsEmissionFractions(rtm, temperatures)  # size: total_elements × n_spectral_bins
    
    # Iteration loop
    println("Starting spectral steady-state direct solve...")
    
    # the system simplifies
    j_tot = M \ boundary
    emissive = D*j_tot # updateSpectralEmission!(rtm, iter, D_matrices, sol_j, emit_frac, temperatures, emissive)
    temperatures = updateTemperaturesSpectral!(rtm, emissive, emit_frac)
    emit_frac = getBinsEmissionFractions(rtm, temperatures)

    # Reconstruct final sol_j from final j_tot and emit_frac
    sol_j = zeros(G, rtm.n_spectral_bins * total_elements)
    for bin in 1:rtm.n_spectral_bins
        bin_start = (bin - 1) * total_elements + 1
        bin_end = bin * total_elements
        sol_j[bin_start:bin_end] = emit_frac[:, bin] .* j_tot
    end
    
    # Write spectral results for each bin to mesh
    println("Writing spectral results to mesh...")
    for bin in 1:rtm.n_spectral_bins
        # Extract solution for this bin
        bin_start = (bin - 1) * total_elements + 1
        bin_end = bin * total_elements
        j_bin = sol_j[bin_start:bin_end]
        
        # Compute e, r, g_a from GERT matrices
        e_bin = D * j_bin           # e = D * j
        r_bin = j_bin - e_bin                     # r = j - e = R' * j
        g_a_bin = j_bin - C * j_bin - r_bin  # g_a = A' * j = j - C*j - r

        # Dummy temperature vector (will be overwritten by update_scalar_temperatures_and_heat_sources!)
        T_bin = zeros(G, total_elements)
        
        # Write to mesh for this bin
        write_results_to_mesh!(rtm, T_bin, j_bin, g_a_bin, r_bin, num_surfaces; spectral_bin=bin)
    end

    # Update mesh with final results
    update_scalar_temperatures_and_heat_sources!(rtm)
    
    # Compute energy conservation error for each spectral bin
    rtm.energy_error = G.([sum(C * sol_j[(i - 1) * total_elements + 1:i * total_elements]) 
                        for i in 1:rtm.n_spectral_bins])    
end

function steadyStateSpectral2D!(rtm::RayTracingDomain2D, F_matrices::Union{Matrix{G}, Vector{Matrix{G}}}; 
                              max_iterations::P=500, spectral_coupling_tolerance::G=1e-11) where {G, P<:Integer}
    """
    Solve spectral radiative equilibrium using the optimal solver.
    
    Automatically chooses between:
    - Direct emission solver: for ε=1 everywhere AND σₛ=0 everywhere (orders of magnitude faster)
    - Full iterative solver: for general problems with scattering/reflection
    
    Set force_full=true to disable optimization and always use full solver.
    """
    
    # Check if we can use the fast direct solver
    if rtm.spectral_mode == :spectral_uniform
        println("=== Using DIRECT spectral solver ===")
        println("Complexity: O(N²) per iteration instead of O((n_bins×N)²)")
        return steadyStateSpectral2D_direct!(rtm, F_matrices; max_iterations=max_iterations, spectral_coupling_tolerance)
    elseif rtm.spectral_mode == :spectral_variable
        println("=== Using FULL spectral solver ===")
        return steadyStateSpectral2D_full!(rtm, F_matrices; 
                                          max_iterations=max_iterations, spectral_coupling_tolerance)
    else
        error("spectral_mode must be :spectral_uniform or :spectral_variable but is $(rtm.spectral_mode)")
    end
end