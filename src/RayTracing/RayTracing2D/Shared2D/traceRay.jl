@inline _binValue(x::Number, bin::Integer) = x
@inline _binValue(x::AbstractVector, bin::Integer) = x[bin]

function traceRay(hmesh::RayTracingDomain2D, p_emit::Point2{G}, dir_emit::Point2{G}, 
                   nudge::G, current_coarse_index::P, spectral_bin::P,
                   local_rng::AbstractRNG) where {G, P<:Integer}
    
    if hmesh.uniform_across_bin[spectral_bin] > -0.1
        # Use fast uniform ray tracing - extinction is same for all bins
        first_face = hmesh.fine_mesh[1][1]
        uniform_beta = _binValue(first_face.kappa_g, spectral_bin) + _binValue(first_face.sigma_s_g, spectral_bin)
        return traceRayUniform(hmesh, p_emit, dir_emit, uniform_beta, nudge, current_coarse_index, local_rng)
    else
        # Use variable extinction ray tracing
        return traceRayVariable(hmesh, p_emit, dir_emit, nudge, current_coarse_index, spectral_bin, local_rng)
    end
end

function traceRayUniform(hmesh::RayTracingDomain2D, p_emit::Point2{G}, dir_emit::Point2{G}, 
                           beta::T, nudge, current_coarse_index::P,
                           local_rng::AbstractRNG) where {G, T, P<:Integer}
    grids  = hmesh.fine_grids_opt::Vector{UniformGrid{G}}
    bboxes = hmesh.fine_bboxes_opt::Vector{Vector{BoundingBox2D{G}}}
    cgrid  = hmesh.coarse_grid_opt::UniformGrid{G}
    cbbox  = hmesh.coarse_bboxes_opt::Vector{BoundingBox2D{G}}
    coarse_mesh = hmesh.coarse_mesh
    all_fine    = hmesh.fine_mesh

    point = p_emit
    direction = dir_emit
    ci::P = current_coarse_index

    S = beta > 0 ? -log(rand(local_rng, G)) / beta : G(Inf)

    @inbounds for _ in 1:10_000
        coarse_face = coarse_mesh[ci]
        u_real, u_index = distToSurface2D(point, direction, coarse_face)

        if S < u_real
            # Gas interaction
            @fastmath point = point + (S - nudge) * direction
            fine_index = findFace2D(all_fine[ci], point, grids[ci], bboxes[ci])
            if fine_index === 0
                return 0  # Ray escaped
            end
            return (ci, fine_index, zero(P), point)
            
        elseif coarse_face.solidWalls[u_index]
            # Wall interaction
            @fastmath point = point + (u_real - nudge) * direction
            fine_index = findFace2D(all_fine[ci], point, grids[ci], bboxes[ci])
            if fine_index === 0
                return 0  # Ray escaped
            end
            u_real, u_index = distToSurface2D(point, direction, all_fine[ci][fine_index])
            return (ci, fine_index, u_index, point)
            
        else
            # Cross to next coarse face
            @fastmath point = point + (u_real + nudge) * direction
            S -= u_real
            
            next_coarse_index = findFace2D(coarse_mesh, point, cgrid, cbbox)
            if next_coarse_index === 0
                return 0  # Ray escaped
            end
            ci = next_coarse_index
        end
    end
    
    return 0  # Maximum iterations reached
end

function traceRayVariable(hmesh::RayTracingDomain2D, p_emit::Point2{G}, dir_emit::Point2{G}, 
                           nudge, current_coarse_index::P, spectral_bin::Int,
                           local_rng::AbstractRNG) where {G, P<:Integer}
    # concrete types asserted once, faces from the typed meshes (see traceRayUniform)
    grids  = hmesh.fine_grids_opt::Vector{UniformGrid{G}}
    bboxes = hmesh.fine_bboxes_opt::Vector{Vector{BoundingBox2D{G}}}
    cgrid  = hmesh.coarse_grid_opt::UniformGrid{G}
    cbbox  = hmesh.coarse_bboxes_opt::Vector{BoundingBox2D{G}}
    coarse_mesh = hmesh.coarse_mesh
    all_fine    = hmesh.fine_mesh

    point = p_emit
    direction = dir_emit
    ci::P = current_coarse_index

    # Sample target optical depth for interaction
    target_tau = -log(rand(local_rng, G))
    accumulated_tau = zero(G)

    @inbounds for _ in 1:10_000
        coarse_face = coarse_mesh[ci]
        u_real, u_index = distToSurface2D(point, direction, coarse_face)

        # Get local extinction coefficient for this spectral bin
        fine_index = findFace2D(all_fine[ci], point, grids[ci], bboxes[ci])
        if fine_index === 0
            return 0  # Ray escaped
        end
        current_fine_face = all_fine[ci][fine_index]

        # Extract extinction for specific spectral bin
        local_beta = _binValue(current_fine_face.kappa_g, spectral_bin) + _binValue(current_fine_face.sigma_s_g, spectral_bin)
        
        # Calculate optical depth to boundary
        tau_to_boundary = local_beta * u_real
        
        if accumulated_tau + tau_to_boundary >= target_tau
            # Gas interaction occurs within this cell
            S = (target_tau - accumulated_tau) / local_beta
            # Gas interaction
            @fastmath point = point + (S - nudge) * direction
            fine_index = findFace2D(all_fine[ci], point, grids[ci], bboxes[ci])
            if fine_index === 0
                return 0  # Ray escaped
            end
            return (ci, fine_index, zero(P), point)
            
        elseif coarse_face.solidWalls[u_index]
            # Wall interaction
            @fastmath point = point + (u_real - nudge) * direction
            fine_index = findFace2D(all_fine[ci], point, grids[ci], bboxes[ci])
            if fine_index === 0
                return 0  # Ray escaped
            end
            u_real, u_index = distToSurface2D(point, direction, all_fine[ci][fine_index])
            return (ci, fine_index, u_index, point)
            
        else
            # Cross to next coarse face
            @fastmath point = point + (u_real + nudge) * direction
            # Update accumulated optical depth
            accumulated_tau += tau_to_boundary
            
            next_coarse_index = findFace2D(coarse_mesh, point, cgrid, cbbox)
            if next_coarse_index === 0
                return 0  # Ray escaped
            end
            ci = next_coarse_index
        end
    end
    
    return 0  # Maximum iterations reached
end

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