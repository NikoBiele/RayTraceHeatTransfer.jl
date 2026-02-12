# Convenience function to update temperatures and heat sources after all spectral bins are computed
function writeTemperaturesHeatSources!(rtm::RayTracingDomain2D, temperatures)

    # spectral
    for ((coarse_idx, fine_idx, wall_index), surface_idx) in rtm.surface_mapping
        sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
        if sub_face.T_in_w[wall_index] > -0.1
            # Compute heat source for prescribed temperature
            sub_face.T_w[wall_index] = sub_face.T_in_w[wall_index]
            sub_face.q_w[wall_index] = sum(sub_face.e_w[wall_index]) - sum(sub_face.g_a_w[wall_index])
        else
            # set computed temperature for prescribed source
            sub_face.T_w[wall_index] = temperatures[surface_idx]
            sub_face.q_w[wall_index] = sub_face.q_in_w[wall_index]
        end
    end
    N_surfs = length(rtm.surface_mapping)
    for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
        sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
        if sub_face.T_in_g > -0.1
            # Compute heat source for prescribed temperature
            sub_face.T_g = sub_face.T_in_g
            sub_face.q_g = sum(sub_face.e_g) - sum(sub_face.g_a_g)
        else
            # set computed temperature for prescribed source
            sub_face.T_g = temperatures[N_surfs + volume_idx]
            sub_face.q_g = sub_face.q_in_g
        end
    end
end

# Convenience function to update temperatures and heat sources after all spectral bins are computed
function writeTemperaturesHeatSources!(domain::ViewFactorDomain3D, temperatures)

    surf_count = 0
    for superface in domain.facesMesh
        for sub_face in superface.subFaces
            surf_count += 1        
            if sub_face.T_in_w > -0.1
                # Compute heat source for prescribed temperature
                sub_face.T_w = sub_face.T_in_w
                sub_face.q_w = sum(sub_face.e_w) - sum(sub_face.g_a_w)
            else
                # set computed temperature for prescribed source
                sub_face.T_w = temperatures[surf_count]
                sub_face.q_w = sub_face.q_in_w
            end
        end
    end
end