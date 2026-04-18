# Updated constructor from existing RayTracingMesh - now with spectral support
function RayTracingDomain2D(rtm::IntermediateMesh2D, verbose::Bool)
    # Extract the type parameters from the input
    VPF = typeof(rtm.coarse_mesh)
    VVPF = typeof(rtm.fine_mesh)
    MT = rtm.F_raw isa Vector ? typeof(rtm.F_raw[1]) : typeof(rtm.F_raw)
    GRID = typeof(rtm.coarse_grid)
    
    # Build optimized cache structures
    
    # 1. Coarse face cache (direct vector access)
    coarse_face_cache = [face for face in rtm.coarse_mesh]
    
    # 2. Fine face cache (flattened structure) 
    fine_face_cache = [Vector{eltype(submesh)}([face for face in submesh]) for submesh in rtm.fine_mesh]
    
    # 3. Pre-compute geometric data for coarse mesh
    coarse_wall_normals = [copy(face.inwardNormals) for face in rtm.coarse_mesh]
    coarse_wall_midpoints = [copy(face.wallMidPoints) for face in rtm.coarse_mesh]
    
    coarse_bounding_boxes = Vector{Tuple{Point2, Point2}}(undef, length(rtm.coarse_mesh))
    for (i, face) in enumerate(rtm.coarse_mesh)
        vertices = face.vertices
        min_x = minimum(v[1] for v in vertices)
        max_x = maximum(v[1] for v in vertices)
        min_y = minimum(v[2] for v in vertices)
        max_y = maximum(v[2] for v in vertices)
        coarse_bounding_boxes[i] = (Point2(min_x, min_y), Point2(max_x, max_y))
    end
    
    # 4. Pre-compute geometric data for fine mesh
    fine_wall_normals = Vector{Vector{Vector{Point2}}}(undef, length(rtm.fine_mesh))
    fine_wall_midpoints = Vector{Vector{Vector{Point2}}}(undef, length(rtm.fine_mesh))
    fine_bounding_boxes = Vector{Vector{Tuple{Point2, Point2}}}(undef, length(rtm.fine_mesh))
    
    for (i, submesh) in enumerate(rtm.fine_mesh)
        fine_wall_normals[i] = [copy(face.inwardNormals) for face in submesh]
        fine_wall_midpoints[i] = [copy(face.wallMidPoints) for face in submesh]
        
        fine_bounding_boxes[i] = Vector{Tuple{Point2, Point2}}(undef, length(submesh))
        for (j, face) in enumerate(submesh)
            vertices = face.vertices
            min_x = minimum(v[1] for v in vertices)
            max_x = maximum(v[1] for v in vertices)
            min_y = minimum(v[2] for v in vertices)
            max_y = maximum(v[2] for v in vertices)
            fine_bounding_boxes[i][j] = (Point2(min_x, min_y), Point2(max_x, max_y))
        end
    end
    
    # 5. Determine spectral mode from the first face
    first_face = fine_face_cache[1][1]
    is_spectral = isa(first_face.kappa_g, Vector)
    n_spectral_bins = is_spectral ? length(first_face.kappa_g) : 1
    
    # 6. Compute mappings
    surface_mapping = Dict{Tuple{Int,Int,Int}, Int}()
    volume_mapping = Dict{Tuple{Int,Int}, Int}()
    surface_index = 1
    volume_index = 1
    surface_areas = []
    volumes = []
    for (coarse_index, coarse_face) in enumerate(rtm.coarse_mesh)
        for (fine_index, fine_face) in enumerate(rtm.fine_mesh[coarse_index])
            for (wall_index, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surface_mapping[(coarse_index, fine_index, wall_index)] = surface_index
                    push!(surface_areas, fine_face.area[wall_index])
                    surface_index += 1
                end
            end
            volume_mapping[(coarse_index, fine_index)] = volume_index
            push!(volumes, fine_face.volume)
            volume_index += 1
        end
    end
    VT = typeof(surface_areas)
    DIII = typeof(surface_mapping)
    DII = typeof(volume_mapping)
        
    rtm_optim = RayTracingDomain2D{VPF,VVPF,MT,VT,DIII,DII,GRID}(
        rtm.coarse_mesh, rtm.fine_mesh, rtm.coarse_grid, rtm.fine_grids,
        rtm.F_raw, rtm.F_smooth,
        surface_areas, volumes, surface_mapping, volume_mapping,
        :not_yet_set, n_spectral_bins, nothing, # NEW spectral fields
        false,
        coarse_face_cache, fine_face_cache,
        coarse_wall_normals, coarse_wall_midpoints,
        fine_wall_normals, fine_wall_midpoints,
        coarse_bounding_boxes, fine_bounding_boxes,
        # Initialize spatial acceleration as nothing (build later)
        nothing,  # coarse_grid_opt
        nothing,  # coarse_bboxes_opt  
        nothing,  # fine_grids_opt
        nothing,  # fine_bboxes_opt
        nothing   # energy_error
    )

    rtm_optim.uniform_extinction = validateExtinctionUniformity!(rtm_optim; verbose=verbose)

    # Update spectral mode based on uniformity
    if validateSpectralUniformity!(rtm_optim; verbose=verbose) && is_spectral
        rtm_optim.spectral_mode = :spectral_uniform
    elseif is_spectral
        rtm_optim.spectral_mode = :spectral_variable
    else
        rtm_optim.spectral_mode = :grey
    end

    return rtm_optim
end

# Updated constructor that builds from scratch (for new meshes) - now with spectral support  
function RayTracingDomain2D(faces::Vector{PolyVolume2D{G}}, Ndiv::Vector{Tuple{P,P}};
                            verbose::Bool=true) where {G, P<:Integer}
    # First create the standard RayTracingMesh
    verbose && println("Building intermediate mesh...")
    standardMesh = IntermediateMesh2D(faces, Ndiv)
    
    # Then convert to optimized version with spectral support
    verbose && println("Optimizing mesh...")
    optimMesh = RayTracingDomain2D(standardMesh, verbose)

    # Build spatial acceleration
    verbose && println("Building spatial acceleration structures...")
    buildSpatialAcceleration!(optimMesh)

    return optimMesh
end