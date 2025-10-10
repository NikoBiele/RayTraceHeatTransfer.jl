function updateTemperaturesSpectral!(rtm::RayTracingMeshOptim, emissive::Vector{G}, emit_frac::Matrix{G}) where {G}
    
    temperatures = zeros(G, length(rtm.surface_mapping) + length(rtm.volume_mapping))

    # Update temperatures for each element
    for ((coarse_idx, fine_idx, wall_idx), surface_idx) in rtm.surface_mapping
        face = rtm.fine_mesh[coarse_idx][fine_idx]
        if face.T_w[wall_idx] < -0.1  # known source flux element
            epsilon_vals = face.epsilon[wall_idx]
            area_val = face.area[wall_idx]
            weighted_epsilon = sum([epsilon_vals[j] * emit_frac[i,j] for j in 1:rtm.n_spectral_bins])
            face.T_w[wall_idx] = (emissive[surface_idx] / (weighted_epsilon * STEFAN_BOLTZMANN * area_val))^(1/4)
            temperatures[surface_idx] = face.T_w[wall_idx]
        end
    end
    for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
        face = rtm.fine_mesh[coarse_idx][fine_idx]
        if face.T_g < -0.1 # known source flux element
            kappa_vals = face.kappa_g
            volume_val = face.volume
            weighted_kappa = sum([kappa_vals[j] * emit_frac[i,j] for j in 1:rtm.n_spectral_bins])
            face.T_g[num_surfaces + volume_idx] = (emissive[num_surfaces + volume_idx] / (4 * STEFAN_BOLTZMANN * weighted_kappa * volume_val))^(1/4)
            temperatures[num_surfaces + volume_idx] = face.T_g[num_surfaces + volume_idx]
        end
    end

    return temperatures
end