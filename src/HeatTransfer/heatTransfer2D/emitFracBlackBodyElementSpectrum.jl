function emitFracBlackBody_element_spectrum(wavelength_bands, temperature, spectral_pos)
    # Skip invalid temperatures
    if !isfinite(temperature) || temperature <= 0.0
        return 0.0
    end
    
    xi = C2/(wavelength_bands[spectral_pos]*temperature)
        
    # Handle extreme cases to prevent overflow/underflow
    F_0_to_lambda_T = 0.0
    if xi > 50.0  # exp(-50) ≈ 2e-22, negligible
        F_0_to_lambda_T = 0.0
    elseif xi < 1e-8  # For very small xi, F approaches 1
        F_0_to_lambda_T = 1.0
    else        
        # Calculate series with overflow protection
        for m = 1:10
            exp_term = exp(-m*xi)
            
            # Early termination if term becomes negligible
            if exp_term < 1e-15
                break
            end
            
            poly_part = xi^3 + 3*xi^2/m + 6*xi/m^2 + 6/m^3
            term = (exp_term/m) * poly_part
            
            # Check for NaN/Inf before adding
            if isfinite(term)
                F_0_to_lambda_T += term
            end
        end

        F_0_to_lambda_T *= 15/pi^4

        # Final safety clamp
        F_0_to_lambda_T = clamp(F_0_to_lambda_T, 0.0, 1.0)
    end

    return F_0_to_lambda_T
end