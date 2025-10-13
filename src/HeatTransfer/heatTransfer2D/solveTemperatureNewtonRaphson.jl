function solve_temperature_newton_raphson(rtm::RayTracingMeshOptim, element_size, measured_powers, absorption_coeffs; 
                                        initial_temp=1000.0, max_iter=10_000, tolerance=1e-12)
    """
    Solve for temperature given measured spectral emission powers using Newton-Raphson
    
    Args:
        measured_powers: Array of measured emission power in each spectral bin
        absorption_coeffs: Array of absorption coefficients ε_i for each bin
        wavelength_bands: Array of wavelength bin boundaries (length = n_bins + 1)
        initial_temp: Initial temperature guess
        max_iter: Maximum iterations
        tolerance: Convergence tolerance for relative temperature change
    """
    
    T = initial_temp
    total_measured_power = sum(measured_powers)
    
    for iter = 1:max_iter
        # Calculate objective function F(T) = Σ P_measured - Σ [ε_i * σT⁴ * f_band_i(T)]
        F = total_measured_power
        dF_dT = 0.0
        
        for i = 1:rtm.n_spectral_bins
            # Get cumulative blackbody fractions at band boundaries
            if i == 1
                F_lower = 0.0
            else
                F_lower = emitFracBlackBody_element_spectrum(rtm.wavelength_band_limits, T, i-1)
            end
            if i == rtm.n_spectral_bins
                F_upper = 1.0
            else
                F_upper = emitFracBlackBody_element_spectrum(rtm.wavelength_band_limits, T, i)   # Upper boundary
            end
                
            # Get derivatives at band boundaries
            if i == 1
                dF_lower_dT = 0.0
            else
                dF_lower_dT = emitFracBlackBody_element_spectrum_derivative(rtm.wavelength_band_limits, T, i-1)
            end
            if i == rtm.n_spectral_bins
                dF_upper_dT = 0.0
            else
                dF_upper_dT = emitFracBlackBody_element_spectrum_derivative(rtm.wavelength_band_limits, T, i)
            end

            # Band fraction and its derivative
            f_band_i = F_upper - F_lower
            df_band_i_dT = dF_upper_dT - dF_lower_dT
            
            # Predicted emission in this bin
            predicted_power_i = f_band_i * absorption_coeffs[i] * element_size * STEFAN_BOLTZMANN * T^4
            
            # Update objective function
            F -= predicted_power_i
            
            # Update derivative: dF/dT = -σ * Σ[ε_i * (4T³*f_band_i + T⁴*df_band_i/dT)]
            dF_dT -= absorption_coeffs[i] * STEFAN_BOLTZMANN * (4*T^3*f_band_i + T^4*df_band_i_dT)
        end
        
        # Newton-Raphson update
        if abs(dF_dT) < 1e-15
            @warn "Derivative too small, Newton-Raphson may not converge"
            break
        end
        
        delta_T = -F / dF_dT
        T_new = T + delta_T
        
        # Ensure temperature stays positive
        T_new = max(T_new, 10.0)  # Minimum reasonable temperature
        
        # Check convergence
        relative_change = abs(delta_T / T)
        if relative_change < tolerance
            println("Converged in $iter iterations to T = $T_new K")
            return T_new
        end
        
        T = T_new
        
        # Debug output
        if iter <= 3 || iter % 5 == 0
            println("Iter $iter: T = $T K, F = $F, relative_change = $relative_change")
        end
    end
    
    @warn "Newton-Raphson did not converge after $max_iter iterations"
    return T
end