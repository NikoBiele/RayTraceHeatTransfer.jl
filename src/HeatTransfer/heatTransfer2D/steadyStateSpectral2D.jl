function steadyStateSpectral2D!(rtm::RayTracingMeshOptim, F_matrices::Union{Matrix{G}, Vector{Matrix{G}}}; 
                              max_iterations::P=500, wavelength_range::Tuple{P,P}=(-7,-3),
                              return_matrices::Bool=false) where {G, P<:Integer}
    """
    Solve spectral radiative equilibrium using the iterative method from your example.
    F_matrices can be either a single matrix (grey/uniform) or vector of matrices (variable spectral)
    """
    
    # Get system matrices
    surface_mapping, volume_mapping = rtm.surface_mapping, rtm.volume_mapping
    num_surfaces = length(surface_mapping)
    num_volumes = length(volume_mapping)
    total_elements = num_surfaces + num_volumes
    
    # Build system matrices based on spectral mode
    D_matrices = Vector{Matrix{G}}()
    C_matrices = Vector{Matrix{G}}()
    M_matrices = Vector{Matrix{G}}()
    if rtm.spectral_mode == :spectral_uniform
        # uniform spectral - single F matrix for all bins
        B, S_infty, A, R, C, D, M = steadyStateGrey2D!(rtm, F_matrices; return_matrices=true, spectral_bin=1)
        # reset the grey energy conservation error (will be re-computed for spectral later)
        rtm.energy_error = nothing

        # For uniform spectral, same matrices apply to all bins
        for bin in 1:rtm.n_spectral_bins
            D_matrices = push!(D_matrices, D)
            C_matrices = push!(C_matrices, C)
            M_matrices = push!(M_matrices, M)
        end
    elseif rtm.spectral_mode == :spectral_variable
        # Variable spectral - use different F matrix for each bin
        
        for bin in 1:rtm.n_spectral_bins
            B, S_infty, A, R, C, D, M = steadyStateGrey2D!(rtm, F_matrices[bin]; return_matrices=true, spectral_bin=bin)
            push!(D_matrices, D)
            push!(C_matrices, C)
            push!(M_matrices, M)
        end

        # reset the grey energy conservation error (will be re-computed for spectral later)
        rtm.energy_error = nothing
    end
    
    # Build block matrix structure
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
    
    # Set boundary conditions from mesh
    boundary, temperatures, emissive = setup_boundary_conditions(rtm, F_matrices, wavelength_range)
    println("Initial temperatures: $temperatures K")
    
    # Iteration loop
    sol_j = zeros(G, rtm.n_spectral_bins * total_elements)
    record_convergence = G[]
    previous_sol_j = zeros(G, rtm.n_spectral_bins * total_elements)
    wavelength_bands = G.(wavelength_band_splits(rtm, wavelength_range))
    emit_frac = getBinsEmissionFractions(rtm, wavelength_bands, temperatures)

    println("Starting spectral steady-state iteration...")
    for iter = 1:max_iterations

        # calculate emissive powers and temperatures
        emissive = updateSpectralEmission!(rtm, iter, D_matrices, sol_j, emit_frac, temperatures, emissive)
        temperatures = updateTemperaturesSpectral!(rtm, emissive, emit_frac)
        emit_frac = getBinsEmissionFractions(rtm, wavelength_bands, temperatures)
        println("Iteration $iter: emissive powers = $emissive")
        
        # calculate inverse of b-e matrix
        emissive_pow_vec = reduce(vcat, [emissive for _ in 1:rtm.n_spectral_bins])
        b_e_matrix = Diagonal(reduce(vcat, [boundary; emissive_pow_vec]))
        sol_j .= block_matrix\(b_e_matrix*[ones(G, total_elements); emit_frac[:]])

        # Check convergence
        convergence_error = maximum(abs.(sol_j - previous_sol_j))
        push!(record_convergence, convergence_error)
        previous_sol_j .= sol_j
        
        println("Iteration $iter: convergence error = $convergence_error")
        
        if iter > 1 && convergence_error < 100_000*eps(Float64)
            println("Converged after $iter iterations")
            Ds_combined = hcat([D_matrices[i] for i in 1:rtm.n_spectral_bins]...)
            emissive = max.(Ds_combined * sol_j, 1e-6)
            emissive = updateSpectralEmission!(rtm, iter, D_matrices, sol_j, emit_frac, temperatures, emissive)
            temperatures = updateTemperaturesSpectral!(rtm, emissive, emit_frac)
            rtm.energy_error = G.([sum(C_matrices[i] * sol_j[(i - 1) * total_elements + 1:i * total_elements]) for i in 1:rtm.n_spectral_bins])
            break
        end
        
        if iter == max_iterations
            println("Warning: Maximum iterations reached. Final error = $convergence_error")
            Ds_combined = hcat([D_matrices[i] for i in 1:rtm.n_spectral_bins]...)
            emissive = max.(Ds_combined * sol_j, 1e-6)
            emissive = updateSpectralEmission!(rtm, iter, D_matrices, sol_j, emit_frac, temperatures, emissive)
            temperatures = updateTemperaturesSpectral!(rtm, emissive, emit_frac)
            rtm.energy_error = G.([sum(C_matrices[i] * sol_j[(i - 1) * total_elements + 1:i * total_elements]) for i in 1:rtm.n_spectral_bins])
        end
    end
    
    # Update mesh with final results
    # update_scalar_temperatures_and_heat_sources!(rtm, wavelength_range; temperatures_in=temperatures)
    
end