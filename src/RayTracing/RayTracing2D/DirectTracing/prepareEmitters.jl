function prepare_emitters(rtm::RayTracingMeshOptim, nudge::G, wavelength_range::Tuple{Int,Int}=(-1,-1), spectral_bin::Int=1) where {G}
    emitters = Emitter[]
    total_energy = G(0.0)

    # Get system matrices
    num_surfaces = length(rtm.surface_mapping)
    num_volumes = length(rtm.volume_mapping)
    total_elements = num_surfaces + num_volumes

    if rtm.spectral_mode == :spectral_uniform || rtm.spectral_mode == :spectral_variable # works!
        temps = zeros(G, total_elements)
        for ((coarse_index, fine_index, wall_index), surface_index) in rtm.surface_mapping
            face = rtm.fine_mesh[coarse_index][fine_index]
            if face.T_in_w[wall_index] > -0.1
                temps[surface_index] = face.T_in_w[wall_index]
            end
        end
        for ((coarse_index, fine_index), volume_index) in rtm.volume_mapping      
            face = rtm.fine_mesh[coarse_index][fine_index]
            if face.T_in_g > -0.1
                temps[num_surfaces + volume_index] = face.T_in_g
            end
        end
        # calculate emissive fractions
        wavelength_bands = G.(wavelength_band_splits(rtm, wavelength_range))
        emit_frac = getBinsEmissionFractions(rtm, wavelength_bands, temps)

        # use emissive fractions to weight absortion to get spectral emissive power
        emissive = zeros(G, total_elements)
        for ((coarse_index, fine_index, wall_index), surface_index) in rtm.surface_mapping
            face = rtm.fine_mesh[coarse_index][fine_index]
            if face.T_in_w[wall_index] > -0.1
                weighted_epsilon = sum([face.epsilon[wall_index][i] * emit_frac[surface_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[surface_index] = emit_frac[surface_index, spectral_bin]*weighted_epsilon*face.area[wall_index]*STEFAN_BOLTZMANN*face.T_in_w[wall_index]^4
                if isfinite(emissive[surface_index])
                    push!(emitters, Emitter(:surface, coarse_index, fine_index, wall_index, emissive[surface_index]))
                    total_energy += emissive[surface_index]
                end
            end
        end
        for ((coarse_index, fine_index), volume_index) in rtm.volume_mapping
            face = rtm.fine_mesh[coarse_index][fine_index]
            if face.T_in_g > -0.1
                volume_count = num_surfaces + volume_index
                weighted_kappa = sum([face.kappa_g[i] * emit_frac[volume_count, i] for i in 1:rtm.n_spectral_bins])
                emissive[volume_count] = emit_frac[volume_count, spectral_bin]*4*STEFAN_BOLTZMANN*weighted_kappa*face.volume*face.T_in_g^4
                if isfinite(emissive[volume_count])
                    push!(emitters, Emitter(:volume, coarse_index, fine_index, 0, emissive[volume_count]))
                    total_energy += emissive[volume_count]
                end
            end
        end
    elseif rtm.spectral_mode == :grey # works!
        # get emissive powers directly
        emissive = zeros(G, total_elements)
        for ((coarse_index, fine_index, wall_index), surface_index) in rtm.surface_mapping
            face = rtm.fine_mesh[coarse_index][fine_index]
            emissive[surface_index] = face.epsilon[wall_index]*face.area[wall_index]*STEFAN_BOLTZMANN*face.T_in_w[wall_index]^4
            if isfinite(emissive[surface_index])
                push!(emitters, Emitter(:surface, coarse_index, fine_index, wall_index, emissive[surface_index]))
                total_energy += emissive[surface_index]
            end
        end
        num_surfaces = length(rtm.surface_mapping)
        for ((coarse_index, fine_index), volume_index) in rtm.volume_mapping
            face = rtm.fine_mesh[coarse_index][fine_index]
            volume_count = num_surfaces + volume_index
            emissive[volume_count] = 4*STEFAN_BOLTZMANN*face.kappa_g*face.volume*face.T_in_g^4
            if isfinite(emissive[volume_count])
                push!(emitters, Emitter(:volume, coarse_index, fine_index, 0, emissive[volume_count]))
                total_energy += emissive[volume_count]
            end
        end
    end
    
    if isnan(total_energy)
        error("No energy detected in mesh")
    end

    # Normalize emitter energies
    for emitter in emitters
        if total_energy > 0.0 && isfinite(emitter.energy)
            emitter.energy /= total_energy
        else
            emitter.energy = 0.0
        end
    end
        
    return emitters, total_energy
end