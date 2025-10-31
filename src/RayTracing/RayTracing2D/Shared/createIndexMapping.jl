function create_index_mapping(rtm::RayTracingDomain2D, rays_total::P) where {P<:Integer}
    surface_index = zero(P)
    volume_index = zero(P)
    surface_mapping = Dict{Tuple{P,P,P}, P}()
    volume_mapping = Dict{Tuple{P,P}, P}()
    
    for (coarse_index, coarse_face) in enumerate(rtm.coarse_mesh)
        for (fine_index, fine_face) in enumerate(rtm.fine_mesh[coarse_index])
            volume_index += one(P)
            volume_mapping[(coarse_index, fine_index)] = volume_index
            for (wall_index, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surface_index += one(P)
                    surface_mapping[(coarse_index, fine_index, wall_index)] = surface_index
                end
            end
        end
    end
    
    return surface_mapping, volume_mapping, surface_index, volume_index
end