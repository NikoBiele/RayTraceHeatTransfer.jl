function traceRay(hmesh::RayTracingDomain2D, p_emit::Point2{G}, dir_emit::Point2{G}, 
                   nudge::G, current_coarse_index::P, spectral_bin::P,
                   local_rng::AbstractRNG) where {G, P<:Integer}
    
    if hmesh.uniform_across_bin[spectral_bin] > -0.1
        # Use fast uniform ray tracing - extinction is same for all bins
        first_face = hmesh.fine_mesh[1][1]
        if isa(first_face.kappa_g, Vector)
            uniform_beta = first_face.kappa_g[spectral_bin] + first_face.sigma_s_g[spectral_bin]
        else
            uniform_beta = first_face.kappa_g + first_face.sigma_s_g
        end
        return traceRayUniform(hmesh, p_emit, dir_emit, uniform_beta, nudge, current_coarse_index, local_rng)
    else
        # Use variable extinction ray tracing
        return traceRayVariable(hmesh, p_emit, dir_emit, nudge, current_coarse_index, spectral_bin, local_rng)
    end
end

# Uniform ray tracing (unchanged - extinction is uniform across all bins)
function traceRayUniform(hmesh::RayTracingDomain2D, p_emit::Point2{G}, dir_emit::Point2{G}, 
                           beta::T, nudge, current_coarse_index::P,
                           local_rng::AbstractRNG) where {G, T, P<:Integer}
    point = p_emit
    direction = dir_emit

    S = beta > 0 ? -log(rand(local_rng, G)) / beta : G(Inf)

    @inbounds for _ in 1:10_000
        coarse_face = hmesh.coarse_face_cache[current_coarse_index]
        u_real, u_index = distToSurface2D(point, direction, coarse_face)

        if S < u_real
            # Gas interaction
            @fastmath point = point + (S - nudge) * direction
            fine_index = findFace2D(hmesh.fine_mesh[current_coarse_index], point, 
                                           hmesh.fine_grids_opt[current_coarse_index],
                                           hmesh.fine_bboxes_opt[current_coarse_index])
            if fine_index === 0
                return 0  # Ray escaped
            end
            return (current_coarse_index, fine_index, zero(P), point)
            
        elseif coarse_face.solidWalls[u_index]
            # Wall interaction
            @fastmath point = point + (u_real - nudge) * direction
            fine_index = findFace2D(hmesh.fine_mesh[current_coarse_index], point,
                                           hmesh.fine_grids_opt[current_coarse_index],
                                           hmesh.fine_bboxes_opt[current_coarse_index])
            if fine_index === 0
                return 0  # Ray escaped
            end
            u_real, u_index = distToSurface2D(point, direction, hmesh.fine_mesh[current_coarse_index][fine_index])
            return (current_coarse_index, fine_index, u_index, point)
            
        else
            # Cross to next coarse face
            @fastmath point = point + (u_real + nudge) * direction
            S -= u_real
            
            next_coarse_index = findFace2D(hmesh.coarse_mesh, point,
                                                   hmesh.coarse_grid_opt, 
                                                   hmesh.coarse_bboxes_opt)
            if next_coarse_index === 0
                return 0  # Ray escaped
            end
            current_coarse_index = next_coarse_index
        end
    end
    
    return 0  # Maximum iterations reached
end

# Updated variable ray tracing with spectral bin support
function traceRayVariable(hmesh::RayTracingDomain2D, p_emit::Point2{G}, dir_emit::Point2{G}, 
                           nudge, current_coarse_index::P, spectral_bin::Int,
                           local_rng::AbstractRNG) where {G, P<:Integer}
    point = p_emit
    direction = dir_emit

    # Sample target optical depth for interaction
    target_tau = -log(rand(local_rng, G))
    accumulated_tau = 0.0

    @inbounds for _ in 1:10_000
        coarse_face = hmesh.coarse_face_cache[current_coarse_index]
        u_real, u_index = distToSurface2D(point, direction, coarse_face)

        # Get local extinction coefficient for this spectral bin
        fine_index = findFace2D(hmesh.fine_mesh[current_coarse_index], point,
                                        hmesh.fine_grids_opt[current_coarse_index],
                                        hmesh.fine_bboxes_opt[current_coarse_index])
        if fine_index === 0
            return 0  # Ray escaped
        end
        current_fine_face = hmesh.fine_mesh[current_coarse_index][fine_index]
        
        # Extract extinction for specific spectral bin
        if isa(current_fine_face.kappa_g, Vector)
            local_beta = current_fine_face.kappa_g[spectral_bin] + current_fine_face.sigma_s_g[spectral_bin]
        else
            local_beta = current_fine_face.kappa_g + current_fine_face.sigma_s_g
        end

        # Calculate optical depth to boundary
        tau_to_boundary = local_beta * u_real

        if accumulated_tau + tau_to_boundary >= target_tau
            # Gas interaction occurs within this cell
            S = (target_tau - accumulated_tau) / local_beta
            # Gas interaction
            @fastmath point = point + (S - nudge) * direction
            fine_index = findFace2D(hmesh.fine_mesh[current_coarse_index], point, 
                                           hmesh.fine_grids_opt[current_coarse_index],
                                           hmesh.fine_bboxes_opt[current_coarse_index])
            if fine_index === 0
                return 0  # Ray escaped
            end
            return (current_coarse_index, fine_index, zero(P), point)
            
        elseif coarse_face.solidWalls[u_index]
            # Wall interaction
            @fastmath point = point + (u_real - nudge) * direction
            fine_index = findFace2D(hmesh.fine_mesh[current_coarse_index], point,
                                           hmesh.fine_grids_opt[current_coarse_index],
                                           hmesh.fine_bboxes_opt[current_coarse_index])
            if fine_index === 0
                return 0  # Ray escaped
            end
            u_real, u_index = distToSurface2D(point, direction, hmesh.fine_mesh[current_coarse_index][fine_index])
            return (current_coarse_index, fine_index, u_index, point)
            
        else
            # Cross to next coarse face
            @fastmath point = point + (u_real + nudge) * direction
            # Update accumulated optical depth
            accumulated_tau += tau_to_boundary

            next_coarse_index = findFace2D(hmesh.coarse_mesh, point,
                                                   hmesh.coarse_grid_opt, 
                                                   hmesh.coarse_bboxes_opt)
            if next_coarse_index === 0
                return 0  # Ray escaped
            end
            current_coarse_index = next_coarse_index
        end
    end
    
    return 0  # Maximum iterations reached
end


#    traceRayPath!(seg_cell, seg_len, hmesh, p_emit, dir_emit, nudge, coarse_index,
#                  volume_mapping, num_surfaces)
#
# Walk a ray fine cell by fine cell with no absorption, appending the global
# volume index and geometric length of every cell crossed to `seg_cell` and
# `seg_len`. Returns `(coarse_index, fine_index, wall_index, point)` at the
# first solid wall hit, or `nothing` if the ray escaped (the caller discards
# the segments appended for such a ray).

function traceRayPath!(seg_cell::Vector{Int32}, seg_len::Vector{Float32},
                       hmesh::RayTracingDomain2D, p_emit::Point2{G}, dir_emit::Point2{G},
                       nudge, coarse0::P,
                       volume_mapping::Dict{Tuple{P,P},P}, num_surfaces::P) where {G, P<:Integer}
    grids  = hmesh.fine_grids_opt::Vector{UniformGrid{G}}
    bboxes = hmesh.fine_bboxes_opt::Vector{Vector{BoundingBox2D{G}}}
    cgrid  = hmesh.coarse_grid_opt::UniformGrid{G}
    cbbox  = hmesh.coarse_bboxes_opt::Vector{BoundingBox2D{G}}
    coarse_mesh = hmesh.coarse_mesh
    all_fine    = hmesh.fine_mesh

    point = p_emit
    direction = dir_emit
    ci::P = coarse0
    fine_mesh = all_fine[ci]
    f0 = findFace2D(fine_mesh, point, grids[ci], bboxes[ci])
    f0 === 0 && return 0
    fi::P = f0

    @inbounds for _ in 1:100_000
        fine_face = fine_mesh[fi]
        u_real, u_index = distToSurface2D(point, direction, fine_face)
        push!(seg_cell, Int32(num_surfaces + volume_mapping[(ci, fi)]))
        push!(seg_len, Float32(u_real))

        if fine_face.solidWalls[u_index]
            @fastmath point = point + (u_real - nudge) * direction
            return (ci, fi, u_index, point)
        end

        @fastmath point = point + (u_real + nudge) * direction
        nf = findFace2D(fine_mesh, point, grids[ci], bboxes[ci])
        if nf === 0
            nc = findFace2D(coarse_mesh, point, cgrid, cbbox)
            nc === 0 && return 0
            ci = nc
            fine_mesh = all_fine[ci]
            nf = findFace2D(fine_mesh, point, grids[ci], bboxes[ci])
            nf === 0 && return 0
        end
        fi = nf
    end
    return 0
end