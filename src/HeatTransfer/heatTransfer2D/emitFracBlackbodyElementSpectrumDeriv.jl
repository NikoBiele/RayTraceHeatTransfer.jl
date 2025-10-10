function emitFracBlackBody_element_spectrum_derivative(wavelength_bands, temperature, spectral_pos)
    # Skip invalid temperatures
    if !isfinite(temperature) || temperature <= 0.0
        return 0.0
    end
    
    xi = C2/(wavelength_bands[spectral_pos]*temperature)
    
    # Handle extreme cases
    dF_dT = 0.0
    if xi > 50.0  # exp(-50) ≈ 2e-22, negligible
        dF_dT = 0.0
    elseif xi < 1e-8  # For very small xi, derivative approaches 0
        dF_dT = 0.0
    else        
        # Calculate derivative series
        for m = 1:10
            exp_term = exp(-m*xi)
            
            # Early termination if term becomes negligible
            if exp_term < 1e-15
                break
            end
            
            # Original polynomial part
            poly_part = xi^3 + 3*xi^2/m + 6*xi/m^2 + 6/m^3
            
            # Derivative of polynomial part w.r.t. xi
            dpoly_dxi = 3*xi^2 + 6*xi/m + 6/m^2
            
            # Apply product rule: d/dxi[exp(-m*xi) * poly] = exp(-m*xi) * (dpoly_dxi - m*poly)
            term_derivative = (exp_term/m) * (dpoly_dxi - m*poly_part)
            
            # Convert dF/dxi to dF/dT using chain rule: dF/dT = (dF/dxi) * (dxi/dT)
            # where dxi/dT = -xi/T
            term_dF_dT = term_derivative * (-xi/temperature)
            
            # Check for NaN/Inf before adding
            if isfinite(term_dF_dT)
                dF_dT += term_dF_dT
            end
        end

        dF_dT *= 15/pi^4
    end

    return dF_dT
end
