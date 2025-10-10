# Replace your find_face_uniform_grid function with this more flexible version:
@inline function find_face_uniform_grid(faces::Vector{PolyFace2D{G}}, point, 
                                        grid) where {G}
    # Convert point to the right type if needed
    pt = Point2{G}(point[1], point[2])
    
    # Find grid cell
    rel_x = pt[1] - grid.origin[1]
    rel_y = pt[2] - grid.origin[2]
    
    cell_i = floor(Int, rel_x * grid.inv_cell_size) + 1
    cell_j = floor(Int, rel_y * grid.inv_cell_size) + 1
    
    # Check bounds
    if cell_i < 1 || cell_i > grid.nx || cell_j < 1 || cell_j > grid.ny
        return nothing
    end
    
    # Test faces in this cell
    @inbounds for face_idx in grid.cells[cell_i, cell_j]
        if point_in_polygon_fast(pt, faces[face_idx])
            return face_idx
        end
    end
    
    return nothing
end

# Replace your find_face_with_bbox_prefilter function with this more flexible version:
@inline function find_face_with_bbox_prefilter(faces::Vector{PolyFace2D{G}}, point, 
                                              bboxes) where {G}
    # Convert point to the right type if needed
    pt = Point2{G}(point[1], point[2])
    
    @inbounds for i in eachindex(faces)
        # Quick bounding box test first
        if point_in_bbox(pt, bboxes[i])
            # Only do expensive point-in-polygon test if bbox test passes
            if point_in_polygon_fast(pt, faces[i])
                return i
            end
        end
    end
    return nothing
end

# Replace your current find_face_optimized function with this more flexible version:
@inline function find_face_optimized(faces::Vector{PolyFace2D{G}}, point, 
                                    grid, bboxes) where {G}
    # Convert point to the right type if needed
    pt = Point2{G}(point[1], point[2])
    
    # Try uniform grid first (fastest for regular grids)
    if grid !== nothing
        result = find_face_uniform_grid(faces, pt, grid)
        if result !== nothing
            return result
        end
    end
    
    # Fallback to bounding box prefilter + linear search
    if bboxes !== nothing
        return find_face_with_bbox_prefilter(faces, pt, bboxes)
    end
    
    # Last resort: return nothing
    return nothing
end

# Fast point-in-bounding-box test
@inline function point_in_bbox(point::Point2{T}, bbox::BoundingBox{T}) where T
    px, py = point[1], point[2]
    return bbox.min_x <= px <= bbox.max_x && bbox.min_y <= py <= bbox.max_y
end

# Optimized point-in-polygon test
@inline function point_in_polygon_fast(point::Point2{G}, face::PolyFace2D{G}) where {G}
    vertices = face.vertices
    n = length(vertices)
    inside = false
    
    px, py = point[1], point[2]
    
    @inbounds begin
        j = n
        for i in 1:n
            xi, yi = vertices[i][1], vertices[i][2]
            xj, yj = vertices[j][1], vertices[j][2]
            
            if ((yi > py) != (yj > py))
                slope = (xj - xi) / (yj - yi)
                intersect_x = xi + slope * (py - yi)
                if px < intersect_x
                    inside = !inside
                end
            end
            j = i
        end
    end
    
    return inside
end