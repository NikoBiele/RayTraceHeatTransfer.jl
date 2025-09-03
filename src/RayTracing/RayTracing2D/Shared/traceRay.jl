function trace_ray(hmesh::RayTracingMeshOptim, p_emit::Point2{G}, dir_emit::Point2{G}, 
                   nudge, current_coarse_index::P) where {G, P<:Integer}
    
    if hmesh.uniform_extinction
        # Use fast uniform ray tracing
        first_face = hmesh.fine_mesh[1][1]  # All faces have same values due to validation
        uniform_beta = first_face.kappa_g + first_face.sigma_s_g
        return trace_ray_uniform(hmesh, p_emit, dir_emit, uniform_beta, nudge, current_coarse_index)
    else
        # Use variable extinction ray tracing
        return trace_ray_variable(hmesh, p_emit, dir_emit, nudge, current_coarse_index)
    end
end

# Uniform ray tracing for performance when all extinction is the same
function trace_ray_uniform(hmesh::RayTracingMeshOptim, p_emit::Point2{G}, dir_emit::Point2{G}, 
                           beta::T, nudge, current_coarse_index::P) where {G, T, P<:Integer}
    point = p_emit
    direction = dir_emit

    S = beta > 0 ? G(-log(rand()) / beta) : G(Inf)

    @inbounds for _ in 1:10_000
        coarse_face = hmesh.coarse_face_cache[current_coarse_index]
        u_real, u_index = distToSurface(point, direction, coarse_face)

        if S < u_real
            # Gas interaction
            @fastmath point = point + (S - nudge) * direction
            fine_index = find_face_optimized(hmesh.fine_mesh[current_coarse_index], point, 
                                           hmesh.fine_grids_opt[current_coarse_index],
                                           hmesh.fine_bboxes_opt[current_coarse_index])
            if fine_index === nothing
                return nothing  # Ray escaped
            end
            return (current_coarse_index, fine_index, zero(P), point)
            
        elseif coarse_face.solidWalls[u_index]
            # Wall interaction
            @fastmath point = point + (u_real - nudge) * direction
            fine_index = find_face_optimized(hmesh.fine_mesh[current_coarse_index], point,
                                           hmesh.fine_grids_opt[current_coarse_index],
                                           hmesh.fine_bboxes_opt[current_coarse_index])
            if fine_index === nothing
                return nothing  # Ray escaped
            end
            u_real, u_index = distToSurface(point, direction, hmesh.fine_mesh[current_coarse_index][fine_index])
            return (current_coarse_index, fine_index, u_index, point)
            
        else
            # Cross to next coarse face
            @fastmath point = point + (u_real + nudge) * direction
            S -= u_real
            
            next_coarse_index = find_face_optimized(hmesh.coarse_mesh, point,
                                                   hmesh.coarse_grid_opt, 
                                                   hmesh.coarse_bboxes_opt)
            if next_coarse_index === nothing
                return nothing  # Ray escaped
            end
            current_coarse_index = next_coarse_index
        end
    end
    
    return nothing  # Maximum iterations reached
end

# Uniform ray tracing for performance when all extinction is the same
function trace_ray_variable(hmesh::RayTracingMeshOptim, p_emit::Point2{G}, dir_emit::Point2{G}, 
                           nudge, current_coarse_index::P) where {G, P<:Integer}
    point = p_emit
    direction = dir_emit

    # Sample target optical depth for interaction
    target_tau = -log(rand())
    accumulated_tau = 0.0

    @inbounds for _ in 1:10_000
        coarse_face = hmesh.coarse_face_cache[current_coarse_index]
        u_real, u_index = distToSurface(point, direction, coarse_face)

        # Get local extinction coefficient directly from face
        fine_index = find_face_optimized(hmesh.fine_mesh[current_coarse_index], point,
                                        hmesh.fine_grids_opt[current_coarse_index],
                                        hmesh.fine_bboxes_opt[current_coarse_index])
        if fine_index === nothing
            return nothing  # Ray escaped
        end
        current_fine_face = hmesh.fine_mesh[current_coarse_index][fine_index]
        local_beta = current_fine_face.kappa_g + current_fine_face.sigma_s_g

        # Calculate optical depth to boundary (only for non-zero extinction)
        tau_to_boundary = local_beta * u_real

        if accumulated_tau + tau_to_boundary >= target_tau
            # Gas interaction occurs within this cell
            S = (target_tau - accumulated_tau) / local_beta
            # Gas interaction
            @fastmath point = point + (S - nudge) * direction
            fine_index = find_face_optimized(hmesh.fine_mesh[current_coarse_index], point, 
                                           hmesh.fine_grids_opt[current_coarse_index],
                                           hmesh.fine_bboxes_opt[current_coarse_index])
            if fine_index === nothing
                return nothing  # Ray escaped
            end
            return (current_coarse_index, fine_index, zero(P), point)
            
        elseif coarse_face.solidWalls[u_index]
            # Wall interaction
            @fastmath point = point + (u_real - nudge) * direction
            fine_index = find_face_optimized(hmesh.fine_mesh[current_coarse_index], point,
                                           hmesh.fine_grids_opt[current_coarse_index],
                                           hmesh.fine_bboxes_opt[current_coarse_index])
            if fine_index === nothing
                return nothing  # Ray escaped
            end
            u_real, u_index = distToSurface(point, direction, hmesh.fine_mesh[current_coarse_index][fine_index])
            return (current_coarse_index, fine_index, u_index, point)
            
        else
            # Cross to next coarse face
            @fastmath point = point + (u_real + nudge) * direction
            # Update accumulated optical depth
            accumulated_tau += tau_to_boundary

            next_coarse_index = find_face_optimized(hmesh.coarse_mesh, point,
                                                   hmesh.coarse_grid_opt, 
                                                   hmesh.coarse_bboxes_opt)
            if next_coarse_index === nothing
                return nothing  # Ray escaped
            end
            current_coarse_index = next_coarse_index
        end
    end
    
    return nothing  # Maximum iterations reached
end