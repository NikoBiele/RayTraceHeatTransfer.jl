function getBinsEmissionFractions(rtm::RayTracingDomain2D, temperatures::Vector{G}) where {G}

    # Calculate spectral emission fractions using Planck function
    emit_frac = zeros(G, length(rtm.surface_mapping)+length(rtm.volume_mapping), rtm.n_spectral_bins)
    num_surfaces = length(rtm.surface_mapping)

    for ((coarse_index, fine_index, wall_index), surface_index) in rtm.surface_mapping
        F_0_to_lambda_T_prev = 0
        for spectral_pos in 1:rtm.n_spectral_bins
            F_0_to_lambda_T = emitFracBlackBody_element_spectrum(rtm.wavelength_band_limits, temperatures[surface_index], spectral_pos)
            if spectral_pos == rtm.n_spectral_bins
                emit_frac[surface_index, spectral_pos] = 1.0 - F_0_to_lambda_T
            else
                emit_frac[surface_index, spectral_pos] = F_0_to_lambda_T - F_0_to_lambda_T_prev
            end
            F_0_to_lambda_T_prev = F_0_to_lambda_T
        end
    end
    for ((coarse_index, fine_index), volume_index) in rtm.volume_mapping
        F_0_to_lambda_T_prev = 0
        for spectral_pos in 1:rtm.n_spectral_bins
            F_0_to_lambda_T = emitFracBlackBody_element_spectrum(rtm.wavelength_band_limits, temperatures[num_surfaces + volume_index], spectral_pos)
            if spectral_pos == rtm.n_spectral_bins
                emit_frac[num_surfaces + volume_index, spectral_pos] = 1.0 - F_0_to_lambda_T
            else
                emit_frac[num_surfaces + volume_index, spectral_pos] = F_0_to_lambda_T - F_0_to_lambda_T_prev
            end
            F_0_to_lambda_T_prev = F_0_to_lambda_T
        end
    end

    return emit_frac
end