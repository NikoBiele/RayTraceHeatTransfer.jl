function setupBoundaryConditions(rtm::RayTracingDomain2D, F_matrices::Union{Matrix{G}, Vector{Matrix{G}}}) where {G}

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
        emitFrac = getBinsEmissionFractions(rtm, temperatures)

        # use emissive fractions to weight absortion to get spectral emissive power
        for ((coarse_face, fine_face, wall_index), surface_index) in rtm.surface_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_w[wall_index] > -0.1
                weighted_epsilon = sum([face.epsilon[wall_index][i] * emitFrac[surface_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[surface_index] = weighted_epsilon*face.area[wall_index]*STEFAN_BOLTZMANN*face.T_in_w[wall_index]^4
                boundary[surface_index] = emissive[surface_index]
            else
                weighted_epsilon = sum([face.epsilon[wall_index][i] * emitFrac[surface_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[surface_index] = weighted_epsilon*face.area[wall_index]*STEFAN_BOLTZMANN*maximum(temperatures)^4
                boundary[surface_index] = face.q_in_w[wall_index]
            end
        end
        for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_g > -0.1
                weighted_kappa = sum([face.kappa_g[i] * emitFrac[num_surfaces + volume_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[num_surfaces + volume_index] = 4*STEFAN_BOLTZMANN*weighted_kappa*face.volume*face.T_in_g^4
                boundary[num_surfaces + volume_index] = emissive[num_surfaces + volume_index]
            else
                weighted_kappa = sum([face.kappa_g[i] * emitFrac[num_surfaces + volume_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[num_surfaces + volume_index] = 4*STEFAN_BOLTZMANN*weighted_kappa*face.volume*maximum(temperatures)^4
                boundary[num_surfaces + volume_index] = face.q_in_g
            end
        end
    elseif rtm.spectral_mode == :spectral_uniform # works!
        
        # calculate emissive fractions
        emitFrac = getBinsEmissionFractions(rtm, temperatures)

        # use emissive fractions to weight absortion to get spectral emissive power
        for ((coarse_face, fine_face, wall_index), surface_index) in rtm.surface_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_w[wall_index] > -0.1
                weighted_epsilon = sum([face.epsilon[wall_index][i] * emitFrac[surface_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[surface_index] = weighted_epsilon*face.area[wall_index]*STEFAN_BOLTZMANN*face.T_in_w[wall_index]^4
                boundary[surface_index] = emissive[surface_index]
            else
                boundary[surface_index] = face.q_in_w[wall_index]
            end
        end
        for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_g > -0.1
                weighted_kappa = sum([face.kappa_g[i] * emitFrac[num_surfaces + volume_index, i] for i in 1:rtm.n_spectral_bins])
                emissive[num_surfaces + volume_index] = 4*STEFAN_BOLTZMANN*weighted_kappa*face.volume*face.T_in_g^4
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

function setupBoundaryConditions(domain::ViewFactorDomain3D{G,P}) where {G,P}
    
    N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    boundary = zeros(G, N_surfs)
    temperatures = zeros(G, N_surfs)
    emissive = zeros(G, N_surfs)
    
    # Get initial temperatures
    surf_count = 0
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            if subface.T_in_w > -0.1
                temperatures[surf_count] = subface.T_in_w
            end
        end
    end
    
    # Calculate emissive fractions
    emitFrac = getBinsEmissionFractions(domain, temperatures)
    
    # Set boundary conditions
    surf_count = 0
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            
            epsilon_vals = isa(subface.epsilon, Vector) ? subface.epsilon : [subface.epsilon]
            weighted_epsilon = sum([epsilon_vals[i] * emitFrac[surf_count, i] 
                                   for i in 1:domain.n_spectral_bins])
            
            if subface.T_in_w > -0.1
                # Known temperature
                emissive[surf_count] = weighted_epsilon * subface.area * 
                                      STEFAN_BOLTZMANN * subface.T_in_w^4
                boundary[surf_count] = emissive[surf_count]
            else
                # Known flux
                emissive[surf_count] = weighted_epsilon * subface.area * 
                                      STEFAN_BOLTZMANN * maximum(temperatures)^4
                boundary[surf_count] = subface.q_in_w
            end
        end
    end
    
    return boundary, temperatures, emissive
end