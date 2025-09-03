function trace_single_ray(hmesh::RayTracingMeshOptim, origin::Point2{Float64},
                        direction::Point2{Float64}, nudge::Float64,
                        current_coarse_index::Int, max_iterations::Int = 100000)
    path = []
    iteration_count = 0
    
    while iteration_count < max_iterations
        iteration_count += 1
        
        # Russian Roulette termination for long rays
        if iteration_count > 1000 && rand() > 0.8
            return nothing
        end
        
        # Check if trace_ray returns nothing (ray escaped)
        ray_result = trace_ray(hmesh, origin, direction, nudge, current_coarse_index)

        if ray_result === nothing
            return nothing
        end
        
        next_coarse_index, next_fine_index, next_wall_index, end_point = ray_result

        if next_wall_index != 0  # Surface interaction
            fine_face = hmesh.fine_mesh[next_coarse_index][next_fine_index]
            epsilon = fine_face.epsilon[next_wall_index]
            
            if rand() < epsilon
                # Check if this wall element is in radiative equilibrium
                if fine_face.T_in_w[next_wall_index] < 0.0
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
            fine_face = hmesh.fine_mesh[next_coarse_index][next_fine_index] ### coarse mesh or fine mesh lookup?

            if rand() < fine_face.sigma_s_g / (fine_face.kappa_g + fine_face.sigma_s_g)
                # Scattering
                direction = isotropicScatter()
                origin = end_point
                push!(path, (next_coarse_index, next_fine_index, -1, :scattering))
            else
                # Check if this gas element is in radiative equilibrium
                fine_face = hmesh.fine_mesh[next_coarse_index][next_fine_index]
                if fine_face.T_in_g < 0.0
                    # Reemission
                    direction = isotropicScatter()
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