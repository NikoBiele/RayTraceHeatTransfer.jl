function setup_boundary_conditions(rtm, F_matrices::Union{Matrix{G}, Vector{Matrix{G}}}, wavelength_range::Tuple{Int,Int}=(-7,-3)) where {G}

    # # Set up surface boundary conditions
    num_surfaces = length(rtm.surface_mapping)
    total_elements = num_surfaces + length(rtm.volume_mapping)
    boundary = zeros(G, total_elements)
    temperatures = zeros(G, total_elements)
    emissive = zeros(G, total_elements)

    for ((coarse_face, fine_face, wall_index), surface_index) in rtm.surface_mapping
        face = rtm.fine_mesh[coarse_face][fine_face]
        if face.T_in_w[wall_index] > -0.1
            temperatures[surface_index] = face.T_in_w[wall_index]
        end
    end
    for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping        
        face = rtm.fine_mesh[coarse_face][fine_face]
        if face.T_in_g > -0.1
            temperatures[num_surfaces + volume_index] = face.T_in_g
        end
    end

    # first get prescribed temperatures and calculate energy
    if rtm.spectral_mode == :spectral_variable # works!
        
        # calculate emissive fractions
        wavelength_bands = G.(wavelength_band_splits(rtm, wavelength_range))
        emit_frac = getBinsEmissionFractions(rtm, wavelength_bands, temperatures)

        # use emissive fractions to weight absortion to get spectral emissive power
        for ((coarse_face, fine_face, wall_index), surface_index) in rtm.surface_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_w[wall_index] > -0.1
                weighted_epsilon = sum([face.epsilon[wall_index][i] * emit_frac[surface_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[surface_index] = weighted_epsilon*face.area[wall_index]*STEFAN_BOLTZMANN*face.T_in_w[wall_index]^4
                boundary[surface_index] = emissive[surface_index]
            else
                weighted_epsilon = sum([face.epsilon[wall_index][i] * emit_frac[surface_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[surface_index] = weighted_epsilon*face.area[wall_index]*STEFAN_BOLTZMANN*maximum(temperatures)^4
                boundary[surface_index] = face.q_in_w[wall_index]
            end
        end
        for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_g > -0.1
                weighted_kappa = sum([face.kappa_g[i] * emit_frac[num_surfaces + volume_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[num_surfaces + volume_index] = 4*STEFAN_BOLTZMANN*weighted_kappa*face.volume*face.T_in_g^4
                boundary[num_surfaces + volume_index] = emissive[num_surfaces + volume_index]
            else
                weighted_kappa = sum([face.kappa_g[i] * emit_frac[num_surfaces + volume_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[num_surfaces + volume_index] = 4*STEFAN_BOLTZMANN*weighted_kappa*face.volume*maximum(temperatures)^4
                boundary[num_surfaces + volume_index] = face.q_in_g
            end
        end
    elseif rtm.spectral_mode == :spectral_uniform # works!
        
        # calculate emissive fractions
        wavelength_bands = G.(wavelength_band_splits(rtm, wavelength_range))
        emit_frac = getBinsEmissionFractions(rtm, wavelength_bands, temperatures)

        # use emissive fractions to weight absortion to get spectral emissive power
        for ((coarse_face, fine_face, wall_index), surface_index) in rtm.surface_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_w[wall_index] > -0.1
                weighted_epsilon = sum([face.epsilon[wall_index][i] * emit_frac[surface_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[surface_index] = weighted_epsilon*face.area[wall_index]*STEFAN_BOLTZMANN*face.T_in_w[wall_index]^4
                boundary[surface_index] = emissive[surface_index]
            else
                boundary[surface_index] = face.q_in_w[wall_index]
            end
        end
        for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_g > -0.1
                weighted_kappa = sum([face.kappa_g[i] * emit_frac[num_surfaces + volume_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[num_surfaces + volume_index] = 4*STEFAN_BOLTZMANN*weighted_kappa*face.Volume*face.T_in_g^4
                boundary[num_surfaces + volume_index] = emissive[num_surfaces + volume_index]
            else
                boundary[num_surfaces + volume_index] = face.q_in_g
            end
        end
    elseif rtm.spectral_mode == :grey # works!
        
        # get emissive powers directly
        for ((coarse_face, fine_face, wall_index), surface_index) in rtm.surface_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_w[wall_index] > -0.1
                emissive[surface_index] = face.epsilon[wall_index]*face.area[wall_index]*STEFAN_BOLTZMANN*face.T_in_w[wall_index]^4
                boundary[surface_index] = emissive[surface_index]
            else
                emissive[surface_index] = STEFAN_BOLTZMANN*maximum(temperatures)^4
                boundary[surface_index] = face.q_in_w[wall_index]
            end        
        end
        for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_g > -0.1
                emissive[num_surfaces + volume_index] = 4*STEFAN_BOLTZMANN*face.kappa_g*face.Volume*face.T_in_g^4
                boundary[num_surfaces + volume_index] = emissive[num_surfaces + volume_index]
            else
                emissive[num_surfaces + volume_index] = STEFAN_BOLTZMANN*maximum(temperatures)^4
                boundary[num_surfaces + volume_index] = face.q_in_g
            end
        end
    end

    return boundary, temperatures, emissive
end