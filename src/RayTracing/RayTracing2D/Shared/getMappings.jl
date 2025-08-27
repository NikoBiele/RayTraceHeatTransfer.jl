function get_mappings(mesh::RayTracingMeshOptim)
    surface_mapping = Dict{Tuple{Int,Int,Int}, Int}()
    volume_mapping = Dict{Tuple{Int,Int}, Int}()
    element_mapping = Dict{Tuple{Int,Int,Int}, Int}()
    
    surface_index = 1
    volume_index = 1

    surface_areas = []
    volumes = []

    for (coarse_index, coarse_face) in enumerate(mesh.coarse_mesh)
        for (fine_index, fine_face) in enumerate(mesh.fine_mesh[coarse_index])
            for (wall_index, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surface_mapping[(coarse_index, fine_index, wall_index)] = surface_index
                    # element_mapping[(coarse_index, fine_index, wall_index)] = surface_index
                    push!(surface_areas, fine_face.area[wall_index])
                    surface_index += 1
                end
            end
            volume_mapping[(coarse_index, fine_index)] = volume_index
            # element_mapping[(coarse_index, fine_index, 0)] = total_surfaces + volume_index
            push!(volumes, fine_face.volume)
            volume_index += 1
        end
    end

    return surface_mapping, volume_mapping, surface_areas, volumes
end