function updateSpectralEmission!(domain::SurfaceDomain3D{G,P}, iter::Int, 
                                   F_matrices::Union{AbstractMatrix, Vector{AbstractMatrix}}, sol_j::Vector{G}, 
                                   emitFrac::AbstractMatrix, temperatures::Vector{G}, 
                                   emissive::Vector{G}) where {G,P}
    
    b = get_b(domain)
    n = size(b,1)
    if iter > 1
        # Update from previous iteration ( using e=sum(D*j) )
        if isa(F_matrices, AbstractVector)
            emissive .= max.(sum([(I-Diagonal(b[:,i])*F_matrices[i]')*sol_j[n*(i-1)+1:i*n] for i in 1:domain.n_spectral_bins]), 10*eps(G))
        elseif isa(F_matrices, AbstractMatrix)
            emissive .= max.(sum([(I-Diagonal(b[:,i])*F_matrices')*sol_j[n*(i-1)+1:i*n] for i in 1:domain.n_spectral_bins]), 10*eps(G))
        else
            error("Incorrect type of 'F_matrices', must be 'Vector{Matrix}' or 'Matrix', got: $(typeof(F_matrices))")
        end
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

function updateSpectralEmission!(rtm::RayTracingDomain2D, iter::Int, F_matrices::Vector{AbstractMatrix},
                                sol_j::Vector{G}, emitFrac::AbstractMatrix, temperatures::Vector{G}, emissive::Vector{G}) where {G}
    num_surfaces = length(rtm.surface_mapping)
    b = get_b(rtm)
    n = size(b,1)

    if iter > 1
        # Update from previous iteration ( using e=sum(D*j) )
        emissive .= max.(sum([(I-Diagonal(b[:,i])*F_matrices[i]')*sol_j[n*(i-1)+1:i*n] for i in 1:rtm.n_spectral_bins]), 10*eps(G))
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