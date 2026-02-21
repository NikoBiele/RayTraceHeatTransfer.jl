function updateSpectralEmission!(domain::ViewFactorDomain3D{G,P}, iter::Int, 
                                   D_matrices::Vector{Matrix{G}}, sol_j::Vector{G}, 
                                   emitFrac::Matrix{G}, temperatures::Vector{G}, 
                                   emissive::Vector{G}) where {G,P}
        
    if iter > 1
        # Update from previous iteration
        Ds_combined = hcat([D_matrices[i] for i in 1:domain.n_spectral_bins]...)
        emissive .= max.(Ds_combined * sol_j, 10*eps(G))
    else
        surf_count = 0
        for superface in domain.facesMesh
            for subface in superface.subFaces
                surf_count += 1
                                
                if subface.T_in_w < -0.1
                    # Known flux - use max temperature
                    emissive[surf_count] = subface.area * STEFAN_BOLTZMANN * maximum(temperatures)^4
                else
                    # Known temperature
                    emissive[surf_count] = subface.area * STEFAN_BOLTZMANN * subface.T_in_w^4
                end
            end
        end
    end
    
    return emissive
end

function updateSpectralEmission!(rtm::RayTracingDomain2D, iter::Int, D_matrices::Vector{Matrix{G}},
                                sol_j::Vector{G}, emitFrac::Matrix{G}, temperatures::Vector{G}, emissive::Vector{G}) where {G}
    num_surfaces = length(rtm.surface_mapping)

    if iter > 1
        # Update temperatures from previous iteration results
        Ds_combined = hcat([D_matrices[i] for i in 1:rtm.n_spectral_bins]...)
        emissive .= max.(Ds_combined * sol_j, 10*eps(G))
    else
        for ((coarse_face, fine_face, wall_index), surface_index) in rtm.surface_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_w[wall_index] < -0.1  # known source flux element
                emissive[surface_index] = face.area[wall_index]*STEFAN_BOLTZMANN*(1+rand())*maximum(temperatures)^4
            else
                emissive[surface_index] = face.area[wall_index]*STEFAN_BOLTZMANN*face.T_in_w[wall_index]^4
            end
        end
        for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if face.T_in_g < -0.1 # known source flux element
                emissive[num_surfaces + volume_index] = 4*STEFAN_BOLTZMANN*face.volume*(1+rand())*maximum(temperatures)^4
            else
                emissive[num_surfaces + volume_index] = 4*STEFAN_BOLTZMANN*face.volume*face.T_in_g^4
            end
        end
    end

    return emissive
    
end