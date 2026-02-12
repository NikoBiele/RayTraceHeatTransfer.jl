# multiple dispatch based on domain type

function getBinsEmissionFractions(rtm::RayTracingDomain2D, temperatures::Vector{G}) where {G}

    # Calculate spectral emission fractions using Planck function
    num_surfaces = length(rtm.surface_mapping)
    emitFrac = zeros(G, num_surfaces+length(rtm.volume_mapping), rtm.n_spectral_bins)

    # for ((coarse_index, fine_index, wall_index), surface_index) in rtm.surface_mapping
    for surface_index in 1:num_surfaces
        F_0_to_lambda_T_prev = 0
        for spectral_pos in 1:rtm.n_spectral_bins
            F_0_to_lambda_T = emitFracBlackBodySpectrum(rtm.wavelength_band_limits, temperatures[surface_index], spectral_pos)
            if spectral_pos == rtm.n_spectral_bins
                emitFrac[surface_index, spectral_pos] = 1.0 - F_0_to_lambda_T
            else
                emitFrac[surface_index, spectral_pos] = F_0_to_lambda_T - F_0_to_lambda_T_prev
            end
            F_0_to_lambda_T_prev = F_0_to_lambda_T
        end
    end
    # for ((coarse_index, fine_index), volume_index) in rtm.volume_mapping
    for volume_index in 1:length(rtm.volume_mapping)
        F_0_to_lambda_T_prev = 0
        for spectral_pos in 1:rtm.n_spectral_bins
            F_0_to_lambda_T = emitFracBlackBodySpectrum(rtm.wavelength_band_limits, temperatures[num_surfaces + volume_index], spectral_pos)
            if spectral_pos == rtm.n_spectral_bins
                emitFrac[num_surfaces + volume_index, spectral_pos] = 1.0 - F_0_to_lambda_T
            else
                emitFrac[num_surfaces + volume_index, spectral_pos] = F_0_to_lambda_T - F_0_to_lambda_T_prev
            end
            F_0_to_lambda_T_prev = F_0_to_lambda_T
        end
    end

    return emitFrac
end

function getBinsEmissionFractions(domain::ViewFactorDomain3D{G,P}, temperatures::Vector{G}) where {G,P}
    
    num_surfaces = sum([length(superface.subFaces) for superface in domain.facesMesh])
    emitFrac = zeros(G, num_surfaces, domain.n_spectral_bins)
    
    for surface_index in 1:num_surfaces
        F_0_to_lambda_T_prev = 0
        for spectral_pos in 1:domain.n_spectral_bins
            F_0_to_lambda_T = emitFracBlackBodySpectrum(domain.wavelength_band_limits, 
                                                                    temperatures[surface_index], 
                                                                    spectral_pos)
            if spectral_pos == domain.n_spectral_bins
                emitFrac[surface_index, spectral_pos] = 1.0 - F_0_to_lambda_T
            else
                emitFrac[surface_index, spectral_pos] = F_0_to_lambda_T - F_0_to_lambda_T_prev
            end
            F_0_to_lambda_T_prev = F_0_to_lambda_T
        end
    end
    
    return emitFrac
end