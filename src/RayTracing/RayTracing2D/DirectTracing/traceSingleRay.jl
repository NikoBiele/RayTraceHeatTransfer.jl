function trace_single_ray(hmesh::RayTracingDomain2D, origin::Point2{G},
                        direction::Point2{G}, nudge::G,
                        current_coarse_index::P, spectral_bin::P=1, max_iterations::P=10_000) where {G, P<:Integer}
    path = []
    iteration_count = 0
    
    while iteration_count < max_iterations
        iteration_count += 1
        
        # Russian Roulette termination for long rays
        if iteration_count > 1000 && rand() > 0.8
            return nothing
        end
        
        # Pass spectral bin to trace_ray
        ray_result = trace_ray(hmesh, origin, direction, nudge, current_coarse_index, spectral_bin)

        if ray_result === nothing
            return nothing
        end
        
        next_coarse_index, next_fine_index, next_wall_index, end_point = ray_result

        if next_wall_index != 0  # Surface interaction
            fine_face = hmesh.fine_mesh[next_coarse_index][next_fine_index]
            epsilon = fine_face.epsilon[next_wall_index]
            
            if rand() < epsilon[spectral_bin]
                # Check if this wall element is in radiative equilibrium for this spectral bin
                wall_temp = fine_face.T_in_w[next_wall_index]
                if wall_temp < 0.0
                    # Reemission
                    p_omit, direction = emit_surface_ray(fine_face, next_wall_index, nudge)
                    # nudge the point a tiny bit towards the midpoint, to ensure we are inside cell
                    origin = end_point + (fine_face.midPoint - end_point) * nudge
                    push!(path, (next_coarse_index, next_fine_index, next_wall_index, :reemission))
                else
                    # True absorption
                    return (:surface, next_coarse_index, next_fine_index, next_wall_index, path)
                end
            else
                # Reflection
                normal = fine_face.outwardNormals[next_wall_index]
                direction = sample_reflection_direction(normal)
                origin = end_point
                push!(path, (next_coarse_index, next_fine_index, next_wall_index, :reflection))
            end
            
        else  # Gas interaction
            fine_face = hmesh.fine_mesh[next_coarse_index][next_fine_index]

            # Get bin-specific scattering properties
            if isa(fine_face.sigma_s_g, Vector)
                local_sigma_s = fine_face.sigma_s_g[spectral_bin]
                local_kappa = fine_face.kappa_g[spectral_bin]
            else
                local_sigma_s = fine_face.sigma_s_g
                local_kappa = fine_face.kappa_g
            end

            if rand() < local_sigma_s / (local_kappa + local_sigma_s)
                # Scattering
                direction = isotropicScatter(nudge) # pass nudge to pass type G
                origin = end_point
                push!(path, (next_coarse_index, next_fine_index, -1, :scattering))
            else
                # Check if this gas element is in radiative equilibrium for this spectral bin
                gas_temp = fine_face.T_in_g
                if gas_temp < 0.0
                    # Reemission
                    direction = isotropicScatter(nudge) # pass nudge to pass type G
                    origin = end_point
                    push!(path, (next_coarse_index, next_fine_index, 0, :reemission))
                else
                    # True absorption
                    return (:gas, next_coarse_index, next_fine_index, 0, path)
                end
            end
        end
        current_coarse_index = next_coarse_index
    end
    
    return nothing
end