function solveTemperatureNewtonRaphson(rtm, element_size, measured_powers, absorption_coeffs;
                                        initial_temp=1000.0, max_iter=10_000, tolerance=1e-12)
    """
    Solve for temperature given measured spectral emission powers using Newton-Raphson.

    Each channel's band fraction is  weight[c] * (F(λ_hi[c]·T) − F(λ_lo[c]·T)), with the first
    coarse bin extended down to λ=0 and the last extended up to λ=∞ —
    identical to the convention in getBinsEmissionFractions, so this inverse
    matches the solver's forward emission model exactly.

    Args:
        measured_powers: Array of measured emission power in each spectral bin/channel
        absorption_coeffs: Array of absorption coefficients ε_i for each bin/channel
        initial_temp: Initial temperature guess
        max_iter: Maximum iterations
        tolerance: Convergence tolerance for relative temperature change

    FIX vs previous version: dF_dT now includes the element_size factor
    (predicted_power_i has it, so the derivative must too; the old Newton step
    was scaled by element_size — same root, perturbed convergence rate).
    """

    T = initial_temp
    total_measured_power = sum(measured_powers)

    for iter = 1:max_iter
        # Calculate objective function F(T) = Σ P_measured - Σ [ε_i * σT⁴ * f_band_i(T)]
        F = total_measured_power
        dF_dT = 0.0

        for i = 1:rtm.n_spectral_bins

            # Cumulative blackbody fractions and derivatives at band boundaries
            if i == 1
                F_lower     = 0.0
                dF_lower_dT = 0.0
            else
                F_lower     = emitFracBlackBodySpectrum(rtm.wavelength_band_limits, T, i-1)
                dF_lower_dT = emitFracBlackBodySpectrumDerivative(rtm.wavelength_band_limits[i-1], T)
            end
            if i == rtm.n_spectral_bins
                F_upper     = 1.0
                dF_upper_dT = 0.0
            else
                F_upper     = emitFracBlackBodySpectrum(rtm.wavelength_band_limits, T, i)
                dF_upper_dT = emitFracBlackBodySpectrumDerivative(rtm.wavelength_band_limits[i], T)
            end
            w = 1.0

            # Band/channel fraction and its derivative
            f_band_i     = w * (F_upper - F_lower)
            df_band_i_dT = w * (dF_upper_dT - dF_lower_dT)

            # Predicted emission in this bin/channel
            predicted_power_i = f_band_i * absorption_coeffs[i] * element_size * STEFAN_BOLTZMANN * T^4

            # Update objective function
            F -= predicted_power_i

            # Update derivative: dF/dT = -σ * Σ[ε_i * A * (4T³*f_band_i + T⁴*df_band_i/dT)]
            dF_dT -= absorption_coeffs[i] * element_size * STEFAN_BOLTZMANN *
                     (4*T^3*f_band_i + T^4*df_band_i_dT)
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
        if iter <= 3 || iter % 20 == 0
            println("Iter $iter: T = $T K, relative_change = $relative_change")
        end
    end

    @warn "Newton-Raphson did not converge after $max_iter iterations"
    return T
end