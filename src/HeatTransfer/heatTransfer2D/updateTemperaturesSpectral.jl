function updateTemperaturesSpectral!(rtm::RayTracingDomain2D, emissive::Vector{G}, emit_frac::Matrix{G}) where {G}
    
    num_surfaces = length(rtm.surface_mapping)
    num_volumes = length(rtm.volume_mapping)
    temperatures = zeros(G, num_surfaces + num_volumes)

    # Update temperatures for surfaces
    for ((coarse_idx, fine_idx, wall_idx), surface_idx) in rtm.surface_mapping
        face = rtm.fine_mesh[coarse_idx][fine_idx]
        
        if face.T_in_w[wall_idx] < -0.1  # Radiative equilibrium (unknown temperature)
            epsilon_vals = face.epsilon[wall_idx]
            area_val = face.area[wall_idx]
            
            # Weighted epsilon based on this element's emission fractions
            weighted_epsilon = sum([epsilon_vals[j] * emit_frac[surface_idx, j] for j in 1:rtm.n_spectral_bins])
            
            # Compute temperature from emissive power
            T_computed = (emissive[surface_idx] / (weighted_epsilon * STEFAN_BOLTZMANN * area_val))^(1/4)
            
            face.T_w[wall_idx] = T_computed
            temperatures[surface_idx] = T_computed
        else
            # Prescribed temperature
            temperatures[surface_idx] = face.T_in_w[wall_idx]
        end
    end
    
    # Update temperatures for volumes
    for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
        face = rtm.fine_mesh[coarse_idx][fine_idx]
        element_idx = num_surfaces + volume_idx
        
        if face.T_in_g < -0.1  # Radiative equilibrium (unknown temperature)
            kappa_vals = face.kappa_g
            volume_val = face.volume
            
            # Weighted kappa based on this element's emission fractions
            weighted_kappa = sum([kappa_vals[j] * emit_frac[element_idx, j] for j in 1:rtm.n_spectral_bins])
            
            # Compute temperature from emissive power
            T_computed = (emissive[element_idx] / (4 * STEFAN_BOLTZMANN * weighted_kappa * volume_val))^(1/4)
            
            face.T_g = T_computed
            temperatures[element_idx] = T_computed
        else
            # Prescribed temperature
            temperatures[element_idx] = face.T_in_g
        end
    end

    return temperatures
end