function trace_ray(hmesh::RayTracingMeshOptim, p_emit::Point2{G}, dir_emit::Point2{G}, 
                   beta::T, nudge, current_coarse_index::P) where {G, T, P<:Integer}
    point = p_emit
    direction = dir_emit

    S = beta > 0 ? G(-log(rand()) / beta) : G(Inf)

    @inbounds for _ in 1:10_000
        coarse_face = hmesh.coarse_face_cache[current_coarse_index]
        u_real, u_index = distToSurface(point, direction, coarse_face)

        if S < u_real
            # Gas interaction - use optimized find_face
            @fastmath point = point + (S - nudge) * direction
            fine_index = find_face_optimized(hmesh.fine_mesh[current_coarse_index], point, 
                                           hmesh.fine_grids_opt[current_coarse_index],
                                           hmesh.fine_bboxes_opt[current_coarse_index])
            if fine_index === nothing
                return nothing  # Ray escaped
            end
            return (current_coarse_index, fine_index, zero(P), point)
            
        elseif coarse_face.solidWalls[u_index]
            # Wall interaction - use optimized find_face
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
            # Cross to next coarse face - use optimized find_face
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