function get_global_index(surface_mapping::Dict{Tuple{P,P,P}, P}, 
                          volume_mapping::Dict{Tuple{P,P}, P}, 
                          num_surfaces::P,
                          coarse_index::P, fine_index::P, wall_index::P, point::Point2{G}) where {G, P<:Integer}
    if wall_index > 0
        index = get(surface_mapping, (coarse_index, fine_index, wall_index), -1)
    else
        index = get(volume_mapping, (coarse_index, fine_index), -1)
        if index != -1
            index += num_surfaces
        end
    end

    return index
end