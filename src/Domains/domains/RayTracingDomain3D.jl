# Constructor from IntermediateMesh3D - with spectral support
function RayTracingDomain3D(domain2d::RayTracingDomain2D, coarse_mesh::Vector{PolyVolume3D{G}}, #fine_mesh::Vector{Vector{Vector{PolyVolume3D{G}}}},
                            F_raw::Union{Matrix{G}, Vector{Matrix{G}}}, F_smooth::Union{Matrix{G}, Vector{Matrix{G}}}) where {G}    
    # 1. Determine spectral mode from the first volume
    first_volume = coarse_mesh[1]
    is_spectral = isa(first_volume.kappa_g, Vector)
    n_spectral_bins = is_spectral ? length(first_volume.kappa_g) : 1
    
    # 2. Compute mappings
    surface_mapping = Dict{Tuple{Int,Int,Int}, Int}()
    # volume_mapping = Dict{Tuple{Int,Int,Int}, Int}()
    surface_index = 1
    # volume_index = 1
    surface_areas = G[]
    # volumes = G[]
    
    for (coarse_index, coarse_volume_2d) in enumerate(coarse_mesh)
        for (coarse_surface_index, coarse_surface_2d) in enumerate(coarse_volume_2d.faces)
            for (fine_surface_index, fine_surface_2d) in enumerate(coarse_surface_2d.subFaces)
                # Map solid faces to global surface indices
                if fine_surface_2d.solidFace !== nothing
                    if fine_surface_2d.solidFace
                        surface_mapping[(coarse_index, coarse_surface_index, fine_surface_index)] = surface_index
                        push!(surface_areas, fine_surface_2d.area)
                        surface_index += 1
                    end
                end
            end
        end
    end
        
    rtm_optim = RayTracingDomain3D{G}(
        coarse_mesh, # fine_mesh,
        F_raw, F_smooth,
        surface_areas, #volumes,
        surface_mapping, # volume_mapping,
        :not_yet_set, n_spectral_bins, nothing,  # Spectral fields
        false,  # uniform_extinction
        # Initialize spatial acceleration as nothing (build later)
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing
    )
    
    # Validate extinction uniformity
    rtm_optim.uniform_extinction = inheritExtinctionUniformity3D!(domain2d) # rtm_optim)
    
    # Validate spectral uniformity
    rtm_optim.spectral_mode = inheritSpectralUniformity3D!(domain2d)
    
    return rtm_optim
end

# Primary user-facing constructor: Extrude 2D domain to 3D
# This creates the full hierarchical 3D mesh from a 2D domain
function RayTracingDomain3D(domain2d::RayTracingDomain2D,
                            z_coords::Vector{G}, Ndim::Int;
                            front_epsilon::Union{Vector{G}, Vector{Vector{G}}},  # One per 2D volume
                            front_q_in::Vector{G},  # One per 2D volume
                            front_T_in::Vector{G},  # One per 2D volume
                            back_epsilon::Union{Vector{G}, Vector{Vector{G}}},   # One per 2D volume
                            back_q_in::Vector{G},   # One per 2D volume
                            back_T_in::Vector{G}    # One per 2D volume
                            ) where {G}
        
    n_coarse = length(domain2d.coarse_mesh)
    # n_z_layers = length(z_coords) - 1
    depth = z_coords[end]
   
    # Step 1: Create coarse 3D mesh (one box per 2D volume, full depth)
    coarse_mesh = Vector{PolyVolume3D{G}}(undef, n_coarse)
    
    for i in 1:n_coarse                
        coarse_mesh[i] = PolyVolume3D{G}(
            domain2d.coarse_mesh[i], depth,
            front_epsilon[i],
            back_epsilon[i],
            front_q_in[i],
            back_q_in[i],
            front_T_in[i],
            back_T_in[i]
        )
    end
    
    # Step 3: Initialize exchange factors (placeholder)
    F_raw = zeros(G, n_coarse, n_coarse)
    F_smooth = zeros(G, n_coarse, n_coarse)

    # Step 4: mesh the coarse faces
    mesh3D = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, n_coarse)
    for i in 1:n_coarse
        mesh3D[i] = Vector{Vector{Vector{Vector{Float64}}}}[]
        for face in coarse_mesh[i].faces
            if length(face.vertices) == 4
                points = [
                    face.vertices[1][1] face.vertices[1][2] face.vertices[1][3];
                    face.vertices[2][1] face.vertices[2][2] face.vertices[2][3];
                    face.vertices[3][1] face.vertices[3][2] face.vertices[3][3];
                    face.vertices[4][1] face.vertices[4][2] face.vertices[4][3]
                ]
                vertex_indices = [1 2 3 4;]
            elseif length(face.vertices) == 3
                points = [
                    face.vertices[1][1] face.vertices[1][2] face.vertices[1][3];
                    face.vertices[2][1] face.vertices[2][2] face.vertices[2][3];
                    face.vertices[3][1] face.vertices[3][2] face.vertices[3][3]
                ]
                vertex_indices = [1 2 3;]
            else
                error("Faces must be triangles or quads")
            end
            face_mesh = meshFaces(points, vertex_indices, Ndim)[1]
            push!(mesh3D[i], face_mesh)
        end
    end
    # num_faces = sum([length(mesh3D[i]) for i in 1:n_coarse])
    # num_points = length(mesh3D[1][1])
    
    # first_epsilon = nothing
    # uniform_epsilon = true

    for i in 1:n_coarse
        for (j, superface) in enumerate(coarse_mesh[i].faces) # previously i
            superface.subFaces = PolyFace3D[]
            
            # First pass: create all subfaces with q=0
            for k in 1:length(mesh3D[i][j][1])
                
                p1 = Point3{G}(mesh3D[i][j][1][k]...)
                p2 = Point3{G}(mesh3D[i][j][2][k]...)
                p3 = Point3{G}(mesh3D[i][j][3][k]...)
                p4 = Point3{G}(mesh3D[i][j][4][k]...)
                
                if isapprox(p3, p4, atol=1e-5) 
                    push!(superface.subFaces, 
                        PolyFace3D([p1, p2, p3],
                        coarse_mesh[i].faces[j].solidFace, 
                        coarse_mesh[i].faces[j].midPoint,
                        coarse_mesh[i].faces[j].epsilon,
                        0.0,
                        coarse_mesh[i].faces[j].T_in_w))
                else
                    push!(superface.subFaces,
                        PolyFace3D([p1, p2, p3, p4],
                        coarse_mesh[i].faces[j].solidFace, 
                        coarse_mesh[i].faces[j].midPoint,
                        coarse_mesh[i].faces[j].epsilon,
                        0.0,
                        coarse_mesh[i].faces[j].T_in_w))
                end
            end
            
            # Second pass: distribute flux proportional to area
            total_area = sum(subface.area for subface in superface.subFaces)
            
            for subface in superface.subFaces
                # Each subface gets flux proportional to its area fraction
                subface.q_in_w = superface.q_in_w * (subface.area / total_area)
                # if first_epsilon === nothing && isa(subface.epsilon, Vector)
                #     first_epsilon = subface.epsilon[1]
                # elseif first_epsilon === nothing && !isa(subface.epsilon, Vector)
                #     first_epsilon = subface.epsilon
                # else
                #     for bin in 1:n_bins
                #         if !isapprox(subface.epsilon[bin], first_epsilon, atol=1e-5)
                #             uniform_epsilon = false
                #         end
                #     end
                # end
            end
        end
    end

    # Step 5: Create mesh
    domain3d = RayTracingDomain3D(domain2d, coarse_mesh, F_raw, F_smooth) 
    
    # reuse 2d acceleration structures, just transfer them (faster than 3d search)
    inheritSpatialAcceleration3D!(domain2d, domain3d)
    
    return domain3d
end

# Reuse 2d spatial acceleration structures
function inheritSpatialAcceleration3D!(domain2d::RayTracingDomain2D, domain3d::RayTracingDomain3D)
    # Optimized cache structures (existing)
    domain3d.from2d_coarse_face_cache = domain2d.coarse_face_cache  # Flattened for direct indexing
    domain3d.from2d_fine_face_cache = domain2d.fine_face_cache  # Pre-allocated fine faces
    
    # Pre-computed geometric data for faster access
    domain3d.from2d_coarse_wall_normals = domain2d.coarse_wall_normals  # Outward normals per face
    domain3d.from2d_coarse_wall_midpoints = domain2d.coarse_wall_midpoints  # Wall midpoints per face
    domain3d.from2d_fine_wall_normals = domain2d.fine_wall_normals  # Fine mesh normals
    domain3d.from2d_fine_wall_midpoints = domain2d.fine_wall_midpoints  # Fine mesh wall midpoints
    
    # Spatial acceleration structures (existing)
    domain3d.from2d_coarse_bounding_boxes = domain2d.coarse_bounding_boxes  # (min, max) per coarse face
    domain3d.from2d_fine_bounding_boxes = domain2d.fine_bounding_boxes  # Bounding boxes for fine faces
    
    # Optimized spatial acceleration structures
    domain3d.from2d_coarse_grid_opt = domain2d.coarse_grid_opt
    domain3d.from2d_coarse_bboxes_opt = domain2d.coarse_bboxes_opt
    domain3d.from2d_fine_grids_opt = domain2d.fine_grids_opt
    domain3d.from2d_fine_bboxes_opt = domain2d.fine_bboxes_opt
    
    return nothing
end