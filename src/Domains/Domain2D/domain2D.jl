const EPSILON = 1e4*eps(Float64)  # Small value for floating-point comparisons

# Define a custom PolyFace2D type
mutable struct PolyFace2D{G<:AbstractFloat, T}

    # geometric variables, fixed types, no uncertainty

    vertices::Vector{Point2{G}}
    solidWalls::Vector{Bool}
    midPoint::Point2{G}
    wallMidPoints::Vector{Point2{G}}
    outwardNormals::Vector{Point2{G}}
    volume::G
    area::Vector{G}
    subFaces::Vector{PolyFace2D{G,T}} # not fixed

    # variables with possibility of uncertainty

    # boundary properties
    epsilon::Vector{T}

    # local gas extinction properties
    kappa_g::T # absorption coefficient [m^-1]
    sigma_s_g::T # scattering coefficient [m^-1]

    # state variables (volume)
    j_g::T # outgoing power [W]
    g_a_g::T # incident absorbed power [W]
    e_g::T # emissive power [W]
    r_g::T # reflected power [W]
    g_g::T # incident power [W]
    i_g::T # total intensity [W*m^(-2)*sr^(-1)]
    q_in_g::T # input source terms [W]
    q_g::T # source terms [W]
    T_in_g::T # input temperatures [K]
    T_g::T # temperatures [K]

    # state variables (walls)
    j_w::Vector{T} # vector of outgoing power [W]
    g_a_w::Vector{T} # vector of incident absorbed power [W]
    e_w::Vector{T} # vector of emissive power [W]
    r_w::Vector{T} # vector of reflected power [W]
    g_w::Vector{T} # vector of incident power [W]
    i_w::Vector{T} # vector of total intensity [W*m^(-2)*sr^(-1)]
    q_in_w::Vector{T} # vector of input source terms [W]
    q_w::Vector{T} # vector of source terms [W]
    T_in_w::Vector{T} # vector of input temperatures [K]
    T_w::Vector{T} # vector of temperatures [K]
end

mutable struct RayTracingMesh{VPF,VVPF,MT,MTU,VT,DIII,DII,BOO,FLOA,GRID}
    coarse_mesh::VPF
    fine_mesh::VVPF
    coarse_grid::GRID  # Not Vector{SpatialGrid}, just a single grid for coarse mesh
    fine_grids::Vector{GRID}  # One grid per coarse face's fine mesh
    F_raw::MT
    F_raw_uncertain::MTU
    F_smooth::MT
    F_smooth_uncertain::MTU
    surface_areas::VT
    volumes::VT
    surface_mapping::DIII
    volume_mapping::DII
    reciprocity_satisfied::BOO
    max_reciprocity_error::FLOA
    conservation_satisfied::BOO
    max_conservation_error::FLOA
    uniform_extinction::Bool
end

struct GridCell{G,P}
    face_indices::Vector{P}
    bounds::Tuple{Point2{G}, Point2{G}}  # min and max points
end

struct SpatialGrid{G,P}
    cells::Matrix{GridCell{G,P}}
    cell_size::Point2{G}
    min_point::Point2{G}
    max_point::Point2{G}
    dims::Tuple{P,P}
end

function RayTracingMesh(faces::Vector{PolyFace2D{G,T}}, Ndiv::Vector{Tuple{P,P}}) where {G<:AbstractFloat, T<:AbstractFloat, P<:Integer}
    meshing_mesh = PolyFace2D{G,T}[]
    for (index, face) in enumerate(faces)
        if length(face.vertices) == 3
            if Ndiv[index][1] != Ndiv[index][2]
                error("Number of divisions must be equal for triangles.")
            else
                push!(meshing_mesh, meshTriangle(face, Ndiv[index][1]))
            end
        elseif length(face.vertices) == 4
            push!(meshing_mesh, meshQuad(face, Ndiv[index][1], Ndiv[index][2]))
        else
            error("Only triangles and quadrilaterals are supported.")
        end
    end
    
    coarse_mesh = faces # unmeshed
    if length(meshing_mesh) == 1
        fine_mesh = [meshing_mesh[1].subFaces] # meshed
    else
        fine_mesh = [submesh.subFaces for submesh in meshing_mesh] # [subsubface for submesh in meshing_mesh for subsubface in submesh.subFaces] # meshed
    end
    
    # Build grids for coarse and fine meshes
    coarse_grid = build_spatial_grid(coarse_mesh, one(P))
    fine_grids = [build_spatial_grid(submesh, one(P)) for submesh in fine_mesh]
    
    # Create appropriate zeros based on input type
    z1 = zeros(T, 2, 2)
    z2 = zeros(T, 2, 2) # uncertainty possibility
    
    rtm = RayTracingMesh(
        coarse_mesh, 
        fine_mesh,
        coarse_grid,
        fine_grids,
        z1,
        z2,
        z1,
        z2,
        zeros(T, 2),
        zeros(T, 2),
        Dict{Tuple{P,P,P}, P}(),
        Dict{Tuple{P,P}, P}(),
        false,
        -one(G), # Float32(-1.0),
        false,
        -one(G), # Float32(-1.0)
        false,
    )

    rtm.uniform_extinction = validate_extinction_consistency!(rtm)

    return rtm
end

# Bounding box for fast rejection tests
struct BoundingBox{T}
    min_x::T
    max_x::T
    min_y::T
    max_y::T
end

# Uniform grid for O(1) spatial queries
struct UniformGrid{T}
    cells::Matrix{Vector{Int}}  # Grid cells containing face indices
    cell_size::T               # Size of each grid cell
    origin::Point2{T}          # Bottom-left corner of grid
    inv_cell_size::T          # Precomputed 1/cell_size for fast division
    nx::Int                   # Number of cells in x direction
    ny::Int                   # Number of cells in y direction
end

# Enhanced RayTracingMesh with optimized data structures
mutable struct RayTracingMeshOptim{VPF,VVPF,MT,MTU,VT,DIII,DII,BOO,FLOA,GRID}
    # Original fields (matching your exact RayTracingMesh structure)
    coarse_mesh::VPF
    fine_mesh::VVPF
    coarse_grid::GRID
    fine_grids::Vector{GRID}
    F_raw::MT
    F_raw_uncertain::MTU
    F_smooth::MT
    F_smooth_uncertain::MTU
    surface_areas::VT
    volumes::VT
    surface_mapping::DIII
    volume_mapping::DII
    reciprocity_satisfied::BOO
    max_reciprocity_error::FLOA
    conservation_satisfied::BOO
    max_conservation_error::FLOA
    uniform_extinction::Bool
    
    # Optimized cache structures
    coarse_face_cache::Vector{PolyFace2D}  # Flattened for direct indexing
    fine_face_cache::Vector{Vector{PolyFace2D}}  # Pre-allocated fine faces
    
    # Pre-computed geometric data for faster access
    coarse_wall_normals::Vector{Vector{Point2}}  # Outward normals per face
    coarse_wall_midpoints::Vector{Vector{Point2}}  # Wall midpoints per face
    fine_wall_normals::Vector{Vector{Vector{Point2}}}  # Fine mesh normals
    fine_wall_midpoints::Vector{Vector{Vector{Point2}}}  # Fine mesh wall midpoints
    
    # Spatial acceleration structures (existing)
    coarse_bounding_boxes::Vector{Tuple{Point2, Point2}}  # (min, max) per coarse face
    fine_bounding_boxes::Vector{Vector{Tuple{Point2, Point2}}}  # Bounding boxes for fine faces
    
    # NEW: Optimized spatial acceleration structures
    coarse_grid_opt::Union{Nothing, UniformGrid}
    coarse_bboxes_opt::Union{Nothing, Vector{BoundingBox}}
    fine_grids_opt::Union{Nothing, Vector{UniformGrid}}
    fine_bboxes_opt::Union{Nothing, Vector{Vector{BoundingBox}}}
end

# Constructor from existing RayTracingMesh - simplified to avoid type parameter issues
function RayTracingMeshOptim(rtm::RayTracingMesh)
    # Extract the type parameters from the input
    VPF = typeof(rtm.coarse_mesh)
    VVPF = typeof(rtm.fine_mesh)
    MT = typeof(rtm.F_raw)
    MTU = typeof(rtm.F_raw_uncertain)
    VT = typeof(rtm.surface_areas)
    DIII = typeof(rtm.surface_mapping)
    DII = typeof(rtm.volume_mapping)
    BOO = typeof(rtm.reciprocity_satisfied)
    FLOA = typeof(rtm.max_reciprocity_error)
    GRID = typeof(rtm.coarse_grid)
    
    # Build optimized cache structures
    println("Building optimized cache structures...")
    
    # 1. Coarse face cache (direct vector access)
    coarse_face_cache = [face for face in rtm.coarse_mesh]
    
    # 2. Fine face cache (flattened structure) 
    fine_face_cache = [Vector{eltype(submesh)}([face for face in submesh]) for submesh in rtm.fine_mesh]
    
    # 3. Pre-compute geometric data for coarse mesh
    println("Pre-computing coarse mesh geometry...")
    coarse_wall_normals = [copy(face.outwardNormals) for face in rtm.coarse_mesh]
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
    println("Pre-computing fine mesh geometry...")
    fine_wall_normals = Vector{Vector{Vector{Point2}}}(undef, length(rtm.fine_mesh))
    fine_wall_midpoints = Vector{Vector{Vector{Point2}}}(undef, length(rtm.fine_mesh))
    fine_bounding_boxes = Vector{Vector{Tuple{Point2, Point2}}}(undef, length(rtm.fine_mesh))
    
    for (i, submesh) in enumerate(rtm.fine_mesh)
        fine_wall_normals[i] = [copy(face.outwardNormals) for face in submesh]
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
    
    println("Optimization cache built successfully!")
    
    return RayTracingMeshOptim{VPF,VVPF,MT,MTU,VT,DIII,DII,BOO,FLOA,GRID}(
        rtm.coarse_mesh, rtm.fine_mesh, rtm.coarse_grid, rtm.fine_grids,
        rtm.F_raw, rtm.F_raw_uncertain, rtm.F_smooth, rtm.F_smooth_uncertain,
        rtm.surface_areas, rtm.volumes, rtm.surface_mapping, rtm.volume_mapping,
        rtm.reciprocity_satisfied, rtm.max_reciprocity_error,
        rtm.conservation_satisfied, rtm.max_conservation_error, rtm.uniform_extinction,
        coarse_face_cache, fine_face_cache,
        coarse_wall_normals, coarse_wall_midpoints,
        fine_wall_normals, fine_wall_midpoints,
        coarse_bounding_boxes, fine_bounding_boxes,
        # Initialize spatial acceleration as nothing (build later)
        nothing,  # coarse_grid_opt
        nothing,  # coarse_bboxes_opt  
        nothing,  # fine_grids_opt
        nothing   # fine_bboxes_opt
    )
end

# Alternative constructor that builds from scratch (for new meshes)
function RayTracingMeshOptim(faces::Vector{PolyFace2D{G,T}}, Ndiv::Vector{Tuple{P,P}}) where {G<:AbstractFloat, T<:AbstractFloat, P<:Integer}
    # First create the standard RayTracingMesh
    standard_mesh = RayTracingMesh(faces, Ndiv)
    
    # Then convert to optimized version
    optimMesh = RayTracingMeshOptim(standard_mesh)

    # Build spatial acceleration
    build_spatial_acceleration!(optimMesh)

    return optimMesh
end

# Triangle constructor with default extinction properties
function PolyFace2D{G,T}(p::SVector{3, Point2{G}}, b::SVector{3,Bool}, 
                         kappa_default::T=zero(T), sigma_s_default::T=zero(T)) where {G<:AbstractFloat, T}
    PolyFace2D(
        # geometric variables
        [p[1], p[2], p[3]], # vertices
        [b[1], b[2], b[3]], # solid walls
        convert(Point2{G}, (p[1]+p[2]+p[3])./3), # mid point
        [convert(Point2{G}, (p[1]+p[2])/2), convert(Point2{G}, (p[2]+p[3])/2), convert(Point2{G}, (p[3]+p[1])/2)], # wall mid points
        
        # outward normals
        [calculate_normal(p[1], p[2], (p[1]+p[2]+p[3])./3),
        calculate_normal(p[2], p[3], (p[1]+p[2]+p[3])./3),
        calculate_normal(p[3], p[1], (p[1]+p[2]+p[3])./3)],
        
        # volume of triangle
        convert(G, 0.5)*(p[1][1]*(p[2][2]-p[3][2])+p[2][1]*(p[3][2]-p[1][2])+p[3][1]*(p[1][2]-p[2][2])),
        
        # area of walls
        [convert(G, norm(p[i]-p[mod(i, 3)+1])) for i = 1:3],
        PolyFace2D{G,T}[], # subfaces

        # boundary properties
        zeros(T, 3),

        # extinction properties
        kappa_default, # kappa_g
        sigma_s_default, # sigma_s_g

        # volume properties
        zero(T), # j_g
        zero(T), # g_a_g
        zero(T), # e_g
        zero(T), # r_g
        zero(T), # g_g
        zero(T), # i_g
        zero(T), # q_in_g
        zero(T), # q_g
        zero(T), # T_in_g
        zero(T), # T_g

        # walls properties
        zeros(T, 3), # j_w
        zeros(T, 3), # g_a_w
        zeros(T, 3), # e_w
        zeros(T, 3), # r_w
        zeros(T, 3), # g_w
        zeros(T, 3), # i_w
        zeros(T, 3), # q_in_w
        zeros(T, 3), # q_w
        zeros(T, 3), # T_in_w
        zeros(T, 3) # T_w
    )
end

# Quadrilateral constructor with default extinction properties
function PolyFace2D{G,T}(p::SVector{4,Point2{G}}, b::SVector{4,Bool}, 
                         kappa_default::T=zero(T), sigma_s_default::T=zero(T)) where {G<:AbstractFloat, T}
    PolyFace2D(
        [p[1], p[2], p[3], p[4]], # vertices
        [b[1], b[2], b[3], b[4]], # solid walls
        convert(Point2{G}, (p[1]+p[2]+p[3]+p[4])/4), # mid point
        [convert(Point2{G}, (p[1]+p[2])/2), 
         convert(Point2{G}, (p[2]+p[3])/2), 
         convert(Point2{G}, (p[3]+p[4])/2), 
         convert(Point2{G}, (p[4]+p[1])/2)], # wall mid points
        
        # outward normals
        [calculate_normal(p[1], p[2], convert(Point2{G}, (p[1]+p[2]+p[3]+p[4])/4)),
         calculate_normal(p[2], p[3], convert(Point2{G}, (p[1]+p[2]+p[3]+p[4])/4)),
         calculate_normal(p[3], p[4], convert(Point2{G}, (p[1]+p[2]+p[3]+p[4])/4)),
         calculate_normal(p[4], p[1], convert(Point2{G}, (p[1]+p[2]+p[3]+p[4])/4))],
        
        # volume
        convert(G, 0.5)*(p[1][1]*(p[2][2]-p[3][2])+p[2][1]*(p[3][2]-p[1][2])+p[3][1]*(p[1][2]-p[2][2])) + 
        convert(G, 0.5)*(p[3][1]*(p[4][2]-p[1][2])+p[4][1]*(p[1][2]-p[3][2])+p[1][1]*(p[3][2]-p[4][2])),
        
        # area of walls
        [convert(G, norm(p[i]-p[mod(i, 4)+1])) for i = 1:4],
        
        PolyFace2D{G,T}[], # subfaces

        # boundary properties
        zeros(T, 4),

        # extinction properties
        kappa_default, # kappa_g
        sigma_s_default, # sigma_s_g

        # volume properties
        zero(T), # j_g
        zero(T), # g_a_g
        zero(T), # e_g
        zero(T), # r_g
        zero(T), # g_g
        zero(T), # i_g
        zero(T), # q_in_g
        zero(T), # q_g
        zero(T), # T_in_g
        zero(T), # T_g

        # walls properties
        zeros(T, 4), # j_w
        zeros(T, 4), # g_a_w
        zeros(T, 4), # e_w
        zeros(T, 4), # r_w
        zeros(T, 4), # g_w
        zeros(T, 4), # i_w
        zeros(T, 4), # q_in_w
        zeros(T, 4), # q_w
        zeros(T, 4), # T_in_w
        zeros(T, 4)  # T_w
    )
end

function calculate_normal(p1::Point2{G}, p2::Point2{G}, midpoint::Point2{G}) where G
    edge = p2 - p1
    normal = normalize(Point2{G}(edge[2], -edge[1])) # normal
    
    # Check if normal is pointing outward
    wall_midpoint = (p1 + p2) / 2
    if dot(normal, wall_midpoint - midpoint) < 0
        normal = -normal  # Flip if pointing inward
    end
    
    return normal
end

# Build uniform grid for a set of faces
function build_uniform_grid(faces::Vector{PolyFace2D{G,T}}, cell_size::G) where {G,T}
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
function compute_bounding_boxes_opt(faces::Vector{PolyFace2D{G,T}}) where {G,T}
    return [BoundingBox(
        minimum(v[1] for v in face.vertices),
        maximum(v[1] for v in face.vertices),
        minimum(v[2] for v in face.vertices),
        maximum(v[2] for v in face.vertices)
    ) for face in faces]
end

# Build optimized spatial structure for a mesh
function build_optimized_spatial_structure(faces::Vector{PolyFace2D{G,T}}) where {G,T}
    if isempty(faces)
        return nothing, nothing
    end
    
    # Estimate optimal grid cell size
    total_area = sum(face.volume for face in faces)
    avg_face_size = sqrt(total_area / length(faces))
    cell_size = avg_face_size * G(2)  # Heuristic: 2x average face size
    
    # Build uniform grid
    grid = build_uniform_grid(faces, cell_size)
    
    # Build bounding boxes for fallback
    bboxes = compute_bounding_boxes_opt(faces)
    
    return grid, bboxes
end

# Main function to build all spatial acceleration structures
function build_spatial_acceleration!(mesh::RayTracingMeshOptim)
    println("Building spatial acceleration structures...")
    
    # Build for coarse mesh
    println("  Building coarse mesh acceleration...")
    mesh.coarse_grid_opt, mesh.coarse_bboxes_opt = build_optimized_spatial_structure(mesh.coarse_mesh)
    
    # Build for each fine mesh
    println("  Building fine mesh acceleration...")
    mesh.fine_grids_opt = Vector{Union{Nothing, UniformGrid}}(undef, length(mesh.fine_mesh))
    mesh.fine_bboxes_opt = Vector{Union{Nothing, Vector{BoundingBox}}}(undef, length(mesh.fine_mesh))
    
    for i in eachindex(mesh.fine_mesh)
        mesh.fine_grids_opt[i], mesh.fine_bboxes_opt[i] = build_optimized_spatial_structure(mesh.fine_mesh[i])
    end
    
    println("Spatial acceleration structures built!")
    return mesh
end

function build_spatial_grid(mesh::Vector{PolyFace2D{G,T}}, ncells::P) where {G<:AbstractFloat,T,P<:Integer}
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
    cells = [GridCell(P[], (Point2(min_x + (i-1)*cell_size[1], min_y + (j-1)*cell_size[2]),
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

function validate_extinction_consistency!(mesh::RayTracingMesh; rtol=1e-10)
    """
    Checks if extinction properties are uniform across all faces.
    Returns true if all faces have approximately the same kappa_g and sigma_s_g values,
    false otherwise. Used to set global uniform_extinction flag.
    """
    
    first_beta = nothing
    
    for (coarse_idx, coarse_face) in enumerate(mesh.coarse_mesh)
        for (fine_idx, fine_face) in enumerate(mesh.fine_mesh[coarse_idx])
            
            # Store first values for comparison
            if first_beta === nothing
                first_beta = fine_face.kappa_g + fine_face.sigma_s_g
            else
                # Check if current values are approximately equal to first values
                if !isapprox(fine_face.kappa_g + fine_face.sigma_s_g, first_beta, rtol=rtol)
                    println("Variable extinction detected - using variable ray tracing")
                    return false
                end
            end
        end
    end
    
    println("Uniform extinction detected (beta_g=$(first_beta)) - using fast uniform ray tracing")
    return true
end