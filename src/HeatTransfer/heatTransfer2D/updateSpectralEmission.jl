function updateSpectralEmission!(rtm::RayTracingMeshOptim, iter::Int, D_matrices::Vector{Matrix{G}},
                                sol_j::Vector{G}, emit_frac::Matrix{G}, temperatures::Vector{G}, emissive::Vector{G}) where {G}
    num_surfaces = length(rtm.surface_mapping)

    # nedenstående virker for spectral uniform, men ikke for spectral variable
    # der skal vi bruge sol_j direkte til at opdatere r og e
    # dernæst benytte Newton-Raphson loopet til at beregne temperaturen
    if iter > 1
        # Update temperatures from previous iteration results
        Ds_combined = hcat([D_matrices[i] for i in 1:rtm.n_spectral_bins]...)
        emissive .= max.(Ds_combined * sol_j, 10*eps(G))
    else
        for ((coarse_face, fine_face, wall_index), surface_index) in rtm.surface_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            # weighted_epsilon = sum([face.epsilon[wall_index][i] * emit_frac[surface_index, i] for i in 1:rtm.n_spectral_bins])
            if face.T_in_w[wall_index] < -0.1  # known source flux element
                emissive[surface_index] = face.area[wall_index]*STEFAN_BOLTZMANN*maximum(temperatures)^4 # weighted_epsilon*
            else
                emissive[surface_index] = face.area[wall_index]*STEFAN_BOLTZMANN*face.T_in_w[wall_index]^4 # weighted_epsilon*
            end
        end
        for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            # weighted_kappa = sum([face.kappa_g[i] * emit_frac[num_surfaces + volume_index, i] for i in 1:rtm.n_spectral_bins])
            if face.T_in_g < -0.1 # known source flux element
                emissive[num_surfaces + volume_index] = 4*STEFAN_BOLTZMANN*face.volume*maximum(temperatures)^4 # weighted_kappa*
            else
                emissive[num_surfaces + volume_index] = 4*STEFAN_BOLTZMANN*face.volume*face.T_in_g^4 # weighted_kappa*
            end
        end
    end

    return emissive
    
end