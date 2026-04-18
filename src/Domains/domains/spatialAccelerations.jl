# Build uniform grid for a set of faces
function buildUniformGrid(faces::Vector{PolyVolume2D{G}}, cell_size::G) where {G}
    if isempty(faces)
        return UniformGrid(Matrix{Vector{Int}}(undef, 1, 1), cell_size, Point2(G(0), G(0)), G(1)/cell_size, 1, 1)
    end
    
    # Find bounding box of all faces
    min_x = min_y = G(Inf)
    max_x = max_y = G(-Inf)
    
    for face in faces
        for vertex in face.vertices
            min_x = min(min_x, vertex[1])
            max_x = max(max_x, vertex[1])
            min_y = min(min_y, vertex[2])
            max_y = max(max_y, vertex[2])
        end
    end
    
    # Add padding
    padding = cell_size * G(0.1)
    min_x -= padding
    min_y -= padding
    max_x += padding
    max_y += padding
    
    # Calculate grid dimensions
    nx = max(1, ceil(Int, (max_x - min_x) / cell_size))
    ny = max(1, ceil(Int, (max_y - min_y) / cell_size))
    
    # Initialize empty grid
    cells = Matrix{Vector{Int}}(undef, nx, ny)
    for i in 1:nx, j in 1:ny
        cells[i, j] = Int[]
    end
    
    # Insert faces into grid cells
    for (face_idx, face) in enumerate(faces)
        # Get face bounding box
        face_min_x = minimum(v[1] for v in face.vertices)
        face_max_x = maximum(v[1] for v in face.vertices)
        face_min_y = minimum(v[2] for v in face.vertices)
        face_max_y = maximum(v[2] for v in face.vertices)
        
        # Find grid cells that intersect with face
        start_i = max(1, floor(Int, (face_min_x - min_x) / cell_size) + 1)
        end_i = min(nx, ceil(Int, (face_max_x - min_x) / cell_size))
        start_j = max(1, floor(Int, (face_min_y - min_y) / cell_size) + 1)
        end_j = min(ny, ceil(Int, (face_max_y - min_y) / cell_size))
        
        # Add face to intersecting cells
        for i in start_i:end_i, j in start_j:end_j
            push!(cells[i, j], face_idx)
        end
    end
    
    origin = Point2(min_x, min_y)
    return UniformGrid(cells, cell_size, origin, G(1)/cell_size, nx, ny)
end

# Compute optimized bounding boxes
function computeBoundingBoxesOpt(faces::Vector{PolyVolume2D{G}}) where {G}
    return [BoundingBox2D(
        minimum(v[1] for v in face.vertices),
        maximum(v[1] for v in face.vertices),
        minimum(v[2] for v in face.vertices),
        maximum(v[2] for v in face.vertices)
    ) for face in faces]
end

# Build optimized spatial structure for a mesh
function buildOptimizedSpatialStructure(faces::Vector{PolyVolume2D{G}}) where {G}
    if isempty(faces)
        return nothing, nothing
    end
    
    # Estimate optimal grid cell size
    total_area = sum(face.volume for face in faces)
    avg_face_size = sqrt(total_area / length(faces))
    cell_size = avg_face_size * G(2)  # Heuristic: 2x average face size
    
    # Build uniform grid
    grid = buildUniformGrid(faces, cell_size)
    
    # Build bounding boxes for fallback
    bboxes = computeBoundingBoxesOpt(faces)
    
    return grid, bboxes
end

# Main function to build all spatial acceleration structures
function buildSpatialAcceleration!(domain::RayTracingDomain2D)
    
    # Build for coarse mesh
    domain.coarse_grid_opt, domain.coarse_bboxes_opt = buildOptimizedSpatialStructure(domain.coarse_mesh)
    
    # Build for each fine domain
    domain.fine_grids_opt = Vector{Union{Nothing, UniformGrid}}(undef, length(domain.fine_mesh))
    domain.fine_bboxes_opt = Vector{Union{Nothing, Vector{BoundingBox2D}}}(undef, length(domain.fine_mesh))
    
    for i in eachindex(domain.fine_mesh)
        domain.fine_grids_opt[i], domain.fine_bboxes_opt[i] = buildOptimizedSpatialStructure(domain.fine_mesh[i])
    end
    
    return domain
end

function buildSpatialGrid(mesh::Vector{PolyVolume2D{G}}, ncells::P) where {G,P<:Integer}
    # Get the underlying numeric type from Measurement type
    numeric_type = if G <: Measurement
        # For Measurement{Float64}, this gives Float64
        G.parameters[1]  # or use: eltype(G)
    else
        G
    end

    # Find bounds of mesh
    min_x = min_y = typemax(numeric_type)
    max_x = max_y = typemin(numeric_type)
    
    for face in mesh
        for vertex in face.vertices
            min_x = min(min_x, vertex[1])
            min_y = min(min_y, vertex[2])
            max_x = max(max_x, vertex[1])
            max_y = max(max_y, vertex[2])
        end
    end
    
    # Add small padding
    padding = (max(max_x - min_x, max_y - min_y)) * 0.01
    min_point = Point2(min_x - padding, min_y - padding)
    max_point = Point2(max_x + padding, max_y + padding)
    
    # Calculate cell size with padding included
    width = max_point[1] - min_point[1]
    height = max_point[2] - min_point[2]
    cell_size = Point2(width / ncells, height / ncells)
    
    # Initialize grid
    cells = [GridCell2D(P[], (Point2(min_x + (i-1)*cell_size[1], min_y + (j-1)*cell_size[2]),
                            Point2(min_x + i*cell_size[1], min_y + j*cell_size[2])))
            for i in 1:ncells, j in 1:ncells]
    
    # Add faces to cells
    for (idx, face) in enumerate(mesh)
        # Find cells that this face overlaps
        min_cell_x = max(1, floor(P, (minimum(v[1] for v in face.vertices) - min_x) / cell_size[1]) + 1)
        max_cell_x = min(ncells, ceil(P, (maximum(v[1] for v in face.vertices) - min_x) / cell_size[1]) + 1)
        min_cell_y = max(1, floor(P, (minimum(v[2] for v in face.vertices) - min_y) / cell_size[2]) + 1)
        max_cell_y = min(ncells, ceil(P, (maximum(v[2] for v in face.vertices) - min_y) / cell_size[2]) + 1)
        
        for i in min_cell_x:max_cell_x, j in min_cell_y:max_cell_y
            push!(cells[i,j].face_indices, idx)
        end
    end
    
    total_indices = sum(length(cell.face_indices) for cell in cells)
    
    SpatialGrid(cells, cell_size, min_point, max_point, (ncells, ncells))
end

# # Build uniform grid from volumes
# function buildUniformGrid3D(volumes::Vector{PolyVolume3D{G}}, 
#                             cell_size::G) where {G}
#     # Find domain bounds
#     all_vertices = reduce(vcat, [vol.vertices for vol in volumes])
    
#     min_x = minimum(v[1] for v in all_vertices)
#     max_x = maximum(v[1] for v in all_vertices)
#     min_y = minimum(v[2] for v in all_vertices)
#     max_y = maximum(v[2] for v in all_vertices)
#     min_z = minimum(v[3] for v in all_vertices)
#     max_z = maximum(v[3] for v in all_vertices)
    
#     origin = Point3{G}(min_x, min_y, min_z)
    
#     # Calculate grid dimensions
#     nx = ceil(Int, (max_x - min_x) / cell_size) + 1
#     ny = ceil(Int, (max_y - min_y) / cell_size) + 1
#     nz = ceil(Int, (max_z - min_z) / cell_size) + 1
    
#     # Initialize empty cells
#     cells = Array{Vector{Int}, 3}(undef, nx, ny, nz)
#     for i in 1:nx, j in 1:ny, k in 1:nz
#         cells[i,j,k] = Int[]
#     end
    
#     # Populate cells with volume indices
#     inv_cell_size = 1 / cell_size
#     for (vol_idx, vol) in enumerate(volumes)
#         # Find grid cells that this volume overlaps
#         vol_min_x = minimum(v[1] for v in vol.vertices)
#         vol_max_x = maximum(v[1] for v in vol.vertices)
#         vol_min_y = minimum(v[2] for v in vol.vertices)
#         vol_max_y = maximum(v[2] for v in vol.vertices)
#         vol_min_z = minimum(v[3] for v in vol.vertices)
#         vol_max_z = maximum(v[3] for v in vol.vertices)
        
#         ix_min = max(1, floor(Int, (vol_min_x - min_x) * inv_cell_size) + 1)
#         ix_max = min(nx, ceil(Int, (vol_max_x - min_x) * inv_cell_size) + 1)
#         iy_min = max(1, floor(Int, (vol_min_y - min_y) * inv_cell_size) + 1)
#         iy_max = min(ny, ceil(Int, (vol_max_y - min_y) * inv_cell_size) + 1)
#         iz_min = max(1, floor(Int, (vol_min_z - min_z) * inv_cell_size) + 1)
#         iz_max = min(nz, ceil(Int, (vol_max_z - min_z) * inv_cell_size) + 1)
        
#         for i in ix_min:ix_max, j in iy_min:iy_max, k in iz_min:iz_max
#             push!(cells[i,j,k], vol_idx)
#         end
#     end
    
#     return UniformGrid3D{G}(cells, cell_size, origin, inv_cell_size, nx, ny, nz)
# end

# # Build bounding boxes for volumes
# function buildBoundingBoxes3D(volumes::Vector{PolyVolume3D{G}}) where {G}
#     bboxes = Vector{BoundingBox3D{G}}(undef, length(volumes))
    
#     for (i, vol) in enumerate(volumes)
#         min_x = minimum(v[1] for v in vol.vertices)
#         max_x = maximum(v[1] for v in vol.vertices)
#         min_y = minimum(v[2] for v in vol.vertices)
#         max_y = maximum(v[2] for v in vol.vertices)
#         min_z = minimum(v[3] for v in vol.vertices)
#         max_z = maximum(v[3] for v in vol.vertices)
        
#         bboxes[i] = BoundingBox3D{G}(min_x, max_x, min_y, max_y, min_z, max_z)
#     end
    
#     return bboxes
# end