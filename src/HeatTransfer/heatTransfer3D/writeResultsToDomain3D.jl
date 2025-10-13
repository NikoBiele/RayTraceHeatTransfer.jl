"""
Write results from exchange factor solution back to 3D domain structure.
Handles grey and spectral modes correctly.
"""

# Grey mode - write scalar results
function write_results_to_domain_grey_3D!(domain::Domain3D_faces, T::Vector{P}, j::Vector{P}, 
                                          Abs::Vector{P}, r::Vector{P}) where P
    surf_count = 0
    
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            
            e_surf = max(j[surf_count] - r[surf_count], 0.0)
            
            # Grey mode - scalar values
            subface.j_w = j[surf_count]
            subface.g_a_w = Abs[surf_count]
            subface.e_w = e_surf
            subface.r_w = r[surf_count]
            subface.g_w = Abs[surf_count] + r[surf_count]
            subface.i_w = j[surf_count] / (π * subface.area)
            if subface.T_in_w < -0.1
                # Prescribed flux
                subface.q_w = subface.q_in_w
                subface.T_w = T[surf_count]
            else
                subface.T_w = subface.T_in_w
                subface.q_w = e_surf - Abs[surf_count]
            end
        end
    end
    
    println("Grey results written: $surf_count surfaces")
end

# Spectral mode - write results for a specific spectral bin
function write_results_to_domain_spectral_bin_3D!(domain::Domain3D_faces, T::Vector{P}, j::Vector{P}, 
                                                  Abs::Vector{P}, r::Vector{P}, spectral_bin::Int) where P
    surf_count = 0
    
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            
            e_surf = max(j[surf_count] - r[surf_count], 0.0)
            
            # Check if spectral vectors already initialized
            if isa(subface.j_w, Vector)
                # Update specific bin
                subface.j_w[spectral_bin] = j[surf_count]
                subface.g_a_w[spectral_bin] = Abs[surf_count]
                subface.e_w[spectral_bin] = e_surf
                subface.r_w[spectral_bin] = r[surf_count]
                subface.g_w[spectral_bin] = Abs[surf_count] + r[surf_count]
                subface.i_w[spectral_bin] = j[surf_count] / (π * subface.area)
            else
                error("Domain surface properties are not spectral vectors but spectral_bin=$spectral_bin was specified")
            end
            
            # Temperature and heat source are ALWAYS scalar (will be updated after all bins)
            # Don't update them here - they will be computed from total energy balance
            if spectral_bin == 1
                # Initialize on first bin
                subface.T_w = T[surf_count]  # Placeholder
                subface.q_w = 0.0  # Will be computed from sum of e_w and g_a_w
            end
        end
    end
    
    println("Spectral bin $spectral_bin results written: $surf_count surfaces")
end

# Main interface - dispatches based on spectral mode
function write_results_to_domain_3D!(domain::Domain3D_faces, T::Vector{P}, j::Vector{P}, 
                                    Abs::Vector{P}, r::Vector{P};
                                    spectral_bin::Union{Nothing,Int}=nothing) where P
    
    if domain.spectral_mode == :grey
        # Grey mode
        if spectral_bin !== nothing && spectral_bin != 1
            @warn "spectral_bin=$spectral_bin specified for grey domain, ignoring"
        end
        write_results_to_domain_grey_3D!(domain, T, j, Abs, r)
        
    elseif domain.spectral_mode == :spectral_uniform || domain.spectral_mode == :spectral_variable
        # Spectral mode
        if spectral_bin === nothing
            error("spectral_bin must be specified for spectral domain")
        end
        write_results_to_domain_spectral_bin_3D!(domain, T, j, Abs, r, spectral_bin)
        
    else
        error("Unknown spectral mode: $(domain.spectral_mode)")
    end
end

# Helper function to update scalar temperatures after all spectral bins computed
function update_scalar_temperatures_3D!(domain::Domain3D_faces{G,P}) where {G,P}
    
    if domain.spectral_mode == :grey
        # Already scalar, nothing to do
        return
    end
    
    for superface in domain.facesMesh
        for subface in superface.subFaces
            if subface.T_in_w < -0.1
                # Radiative equilibrium - use Newton-Raphson to get scalar temperature
                measured_powers = isa(subface.e_w, Vector) ? subface.e_w : [subface.e_w]
                epsilon_vals = isa(subface.epsilon, Vector) ? subface.epsilon : [subface.epsilon]
                
                # Get initial temperature guess
                max_temp = maximum([sf.T_in_w for sf in domain.facesMesh for ssf in sf.subFaces if ssf.T_in_w > 0])
                
                # Solve for scalar temperature
                subface.T_w = solve_temperature_newton_raphson_3D(
                    domain, subface.area, measured_powers, epsilon_vals;
                    initial_temp=max_temp, max_iter=10_000, tolerance=1000*eps(G)
                )
                
                # Prescribed heat source (scalar)
                subface.q_w = subface.q_in_w
            else
                # Prescribed temperature (scalar)
                subface.T_w = subface.T_in_w
                
                # Heat source from total energy balance (scalar)
                total_e = sum(subface.e_w)
                total_abs = sum(subface.g_a_w)
                subface.q_w = total_e - total_abs
            end
        end
    end
end

# Newton-Raphson solver for 3D (adapted from 2D version)
function solve_temperature_newton_raphson_3D(domain::Domain3D_faces, element_size, measured_powers, 
                                            absorption_coeffs; 
                                            initial_temp=1000.0, max_iter=10_000, tolerance=1e-12)
    
    T = initial_temp
    total_measured_power = sum(measured_powers)
    n_bins = domain.n_spectral_bins
    
    for iter = 1:max_iter
        F = total_measured_power
        dF_dT = 0.0
        
        for i = 1:n_bins
            # Get band fractions
            if i == 1
                F_lower = 0.0
                dF_lower_dT = 0.0
            else
                F_lower = emitFracBlackBody_element_spectrum(domain.wavelength_band_limits, T, i-1)
                dF_lower_dT = emitFracBlackBody_element_spectrum_derivative(domain.wavelength_band_limits, T, i-1)
            end
            
            if i == n_bins
                F_upper = 1.0
                dF_upper_dT = 0.0
            else
                F_upper = emitFracBlackBody_element_spectrum(domain.wavelength_band_limits, T, i)
                dF_upper_dT = emitFracBlackBody_element_spectrum_derivative(domain.wavelength_band_limits, T, i)
            end
            
            f_band_i = F_upper - F_lower
            df_band_i_dT = dF_upper_dT - dF_lower_dT
            
            predicted_power_i = f_band_i * absorption_coeffs[i] * element_size * STEFAN_BOLTZMANN * T^4
            F -= predicted_power_i
            dF_dT -= absorption_coeffs[i] * element_size * STEFAN_BOLTZMANN * (4*T^3*f_band_i + T^4*df_band_i_dT)
        end
        
        if abs(dF_dT) < 1e-15
            @warn "Derivative too small, Newton-Raphson may not converge"
            break
        end
        
        delta_T = -F / dF_dT
        T_new = max(T + delta_T, 10.0)
        
        if abs(delta_T / T) < tolerance
            return T_new
        end
        
        T = T_new
    end
    
    @warn "Newton-Raphson did not converge after $max_iter iterations"
    return T
end