const EPSILON = 1e4*eps(Float64)  # Small value for floating-point comparisons

# Define a custom PolyFace2D type
mutable struct PolyFace2D{G}

    # geometric variables, fixed types, no uncertainty (unchanged)
    vertices::Vector{Point2{G}}
    solidWalls::Vector{Bool}
    midPoint::Point2{G}
    wallMidPoints::Vector{Point2{G}}
    outwardNormals::Vector{Point2{G}}
    volume::G
    area::Vector{G}
    subFaces::Vector{PolyFace2D{G}} # not fixed

    # boundary properties (unchanged)
    epsilon::Vector{Union{G, Vector{G}}} # spectral emissivity

    # UPDATED: local gas extinction properties - now Union for spectral support
    kappa_g::Union{G, Vector{G}}     # absorption coefficient [m^-1] - scalar (grey) or vector (spectral)
    sigma_s_g::Union{G, Vector{G}}   # scattering coefficient [m^-1] - scalar (grey) or vector (spectral)

    # UPDATED: state variables (volume) - now Union for spectral support
    j_g::Union{G, Vector{G}}         # outgoing power [W]
    g_a_g::Union{G, Vector{G}}       # incident absorbed power [W]
    e_g::Union{G, Vector{G}}         # emissive power [W]
    r_g::Union{G, Vector{G}}         # reflected power [W]
    g_g::Union{G, Vector{G}}         # incident power [W]
    i_g::Union{G, Vector{G}}         # total intensity [W*m^(-2)*sr^(-1)]
    q_in_g::G      # input source terms [W]
    q_g::G         # source terms [W]
    T_in_g::G      # input temperatures [K]
    T_g::G         # temperatures [K]

    # UPDATED: state variables (walls) - Vector of Union for wall × spectral support
    j_w::Vector{Union{G, Vector{G}}}     # vector of outgoing power [W] - each wall can be grey or spectral
    g_a_w::Vector{Union{G, Vector{G}}}   # vector of incident absorbed power [W]
    e_w::Vector{Union{G, Vector{G}}}     # vector of emissive power [W]
    r_w::Vector{Union{G, Vector{G}}}     # vector of reflected power [W]
    g_w::Vector{Union{G, Vector{G}}}     # vector of incident power [W]
    i_w::Vector{Union{G, Vector{G}}}     # vector of total intensity [W*m^(-2)*sr^(-1)]
    q_in_w::Vector{G}  # vector of input source terms [W]
    q_w::Vector{G}     # vector of source terms [W]
    T_in_w::Vector{G}  # vector of input temperatures [K]
    T_w::Vector{G}     # vector of temperatures [K]
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

mutable struct RayTracingMesh{VPF,VVPF,MT,VT,DIII,DII,GRID}
    coarse_mesh::VPF
    fine_mesh::VVPF
    coarse_grid::GRID  
    fine_grids::Vector{GRID}
    
    # UPDATED: Exchange factor matrices - Union for spectral support
    F_raw::Union{MT, Vector{MT}}              # Single matrix (grey) or vector of matrices (spectral)
    # F_raw_uncertain::Union{MTU, Vector{MTU}}  # Single matrix (grey) or vector of matrices (spectral)
    F_smooth::Union{MT, Vector{MT}}           # Single matrix (grey) or vector of matrices (spectral)
    # F_smooth_uncertain::Union{MTU, Vector{MTU}} # Single matrix (grey) or vector of matrices (spectral)
    
    surface_areas::VT
    volumes::VT
    surface_mapping::DIII
    volume_mapping::DII
    uniform_extinction::Bool
    
    # NEW: Spectral metadata
    spectral_mode::Symbol        # :grey, :spectral_uniform, :spectral_variable
    n_spectral_bins::Int        # Number of spectral bins (1 for grey)
end

# Updated constructor
function RayTracingMesh(faces::Vector{PolyFace2D{G}}, Ndiv::Vector{Tuple{P,P}}) where {G, P<:Integer}
    meshing_mesh = PolyFace2D{G}[]
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
        fine_mesh = [submesh.subFaces for submesh in meshing_mesh] # meshed
    end
    
    # Build grids for coarse and fine meshes
    coarse_grid = build_spatial_grid(coarse_mesh, one(P))
    fine_grids = [build_spatial_grid(submesh, one(P)) for submesh in fine_mesh]
    
    # Determine spectral mode from the first face
    first_face = fine_mesh[1][1]
    is_spectral = isa(first_face.kappa_g, Vector)
    n_bins = is_spectral ? length(first_face.kappa_g) : 1
    
    # Initialize F matrices based on spectral mode
    if is_spectral
        # Spectral mode - create vector of matrices (will be populated during ray tracing)
        F_raw = Matrix{G}[]
        push!(F_raw, zeros(G, 2, 2))
        # F_raw_uncertain = Matrix{G}[]
        # push!(F_raw_uncertain, zeros(G, 2, 2))
        F_smooth = Matrix{G}[]
        push!(F_smooth, zeros(G, 2, 2))
        # F_smooth_uncertain = Matrix{G}[]
        # push!(F_smooth_uncertain, zeros(G, 2, 2))
        spectral_mode = :spectral_variable  # Will be determined by validate_extinction_consistency!
    else
        # Grey mode - single matrices
        F_raw = zeros(G, 2, 2)
        # F_raw_uncertain = zeros(G, 2, 2)
        F_smooth = zeros(G, 2, 2)
        # F_smooth_uncertain = zeros(G, 2, 2)
        spectral_mode = :grey
    end
    
    rtm = RayTracingMesh(
        coarse_mesh, 
        fine_mesh,
        coarse_grid,
        fine_grids,
        F_raw,
        # F_raw_uncertain,
        F_smooth,
        # F_smooth_uncertain,
        Vector{G}(),
        Vector{G}(),
        Dict{Tuple{Int,Int,Int}, Int}(),
        Dict{Tuple{Int,Int}, Int}(),
        false,
        spectral_mode,
        n_bins
    )

    rtm.uniform_extinction = validate_extinction_consistency!(rtm)
    
    # Update spectral mode based on uniformity
    if rtm.uniform_extinction && is_spectral
        rtm.spectral_mode = :spectral_uniform
    end

    return rtm
end

function get_F_matrix(mesh::RayTracingMesh, spectral_bin::Int=1)
    if mesh.spectral_mode == :spectral_variable
        return mesh.F_smooth[spectral_bin]
    else
        # For :grey and :spectral_uniform, use the single matrix
        return mesh.F_smooth
    end
end

function get_F_matrices(mesh::RayTracingMesh)
    if mesh.spectral_mode == :spectral_variable
        return mesh.F_smooth
    else
        # For :grey and :spectral_uniform, return single matrix as vector for consistent interface
        return [mesh.F_smooth]
    end
end

# Bounding box for fast rejection tests
struct BoundingBox{G}
    min_x::G
    max_x::G
    min_y::G
    max_y::G
end

# Uniform grid for O(1) spatial queries
struct UniformGrid{G}
    cells::Matrix{Vector{Int}}  # Grid cells containing face indices
    cell_size::G               # Size of each grid cell
    origin::Point2{G}          # Bottom-left corner of grid
    inv_cell_size::G          # Precomputed 1/cell_size for fast division
    nx::Int                   # Number of cells in x direction
    ny::Int                   # Number of cells in y direction
end

# Enhanced RayTracingMeshOptim with spectral support
mutable struct RayTracingMeshOptim{VPF,VVPF,MT,VT,DIII,DII,GRID}
    # Original fields (matching your exact RayTracingMesh structure)
    coarse_mesh::VPF
    fine_mesh::VVPF
    coarse_grid::GRID
    fine_grids::Vector{GRID}
    F_raw::Union{MT, Vector{MT}}              # UPDATED: Union for spectral support
    # F_raw_uncertain::Union{MTU, Vector{MTU}}  # UPDATED: Union for spectral support
    F_smooth::Union{MT, Vector{MT}}           # UPDATED: Union for spectral support
    # F_smooth_uncertain::Union{MTU, Vector{MTU}} # UPDATED: Union for spectral support
    surface_areas::VT
    volumes::VT
    surface_mapping::DIII
    volume_mapping::DII
    uniform_extinction::Bool
    
    # NEW: Spectral metadata
    spectral_mode::Symbol        # :grey, :spectral_uniform, :spectral_variable
    n_spectral_bins::Int        # Number of spectral bins (1 for grey)
    
    # Optimized cache structures (existing)
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
    
    # Optimized spatial acceleration structures
    coarse_grid_opt::Union{Nothing, UniformGrid}
    coarse_bboxes_opt::Union{Nothing, Vector{BoundingBox}}
    fine_grids_opt::Union{Nothing, Vector{UniformGrid}}
    fine_bboxes_opt::Union{Nothing, Vector{Vector{BoundingBox}}}

    # automatic calculation of global energy conservation error
    energy_error::Union{Nothing, G, Vector{G}} where {G}
end

# Updated constructor from existing RayTracingMesh - now with spectral support
function RayTracingMeshOptim(rtm::RayTracingMesh)
    # Extract the type parameters from the input
    VPF = typeof(rtm.coarse_mesh)
    VVPF = typeof(rtm.fine_mesh)
    MT = rtm.F_raw isa Vector ? typeof(rtm.F_raw[1]) : typeof(rtm.F_raw)
    # MTU = rtm.F_raw_uncertain isa Vector ? typeof(rtm.F_raw_uncertain[1]) : typeof(rtm.F_raw_uncertain)
    VT = typeof(rtm.surface_areas)
    DIII = typeof(rtm.surface_mapping)
    DII = typeof(rtm.volume_mapping)
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
    
    # 5. Determine spectral mode from the first face
    first_face = fine_face_cache[1][1]
    is_spectral = isa(first_face.kappa_g, Vector)
    n_bins = is_spectral ? length(first_face.kappa_g) : 1
    
    # Initialize F matrices based on spectral mode
    if is_spectral
        # Check if extinction is uniform across all bins and faces
        spectral_mode = rtm.uniform_extinction ? :spectral_uniform : :spectral_variable
    else
        spectral_mode = :grey
    end
    
    # 6. Compute mappings
    surface_mapping = Dict{Tuple{Int,Int,Int}, Int}()
    volume_mapping = Dict{Tuple{Int,Int}, Int}()
    element_mapping = Dict{Tuple{Int,Int,Int}, Int}()
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
    
    println("Optimization cache built successfully!")
    
    return RayTracingMeshOptim{VPF,VVPF,MT,VT,DIII,DII,GRID}(
        rtm.coarse_mesh, rtm.fine_mesh, rtm.coarse_grid, rtm.fine_grids,
        rtm.F_raw, # rtm.F_raw_uncertain,
        rtm.F_smooth, # rtm.F_smooth_uncertain,
        surface_areas, volumes, surface_mapping, volume_mapping,
        rtm.uniform_extinction,
        spectral_mode, n_bins,  # NEW spectral fields
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
end

# Updated constructor that builds from scratch (for new meshes) - now with spectral support  
function RayTracingMeshOptim(faces::Vector{PolyFace2D{G}}, Ndiv::Vector{Tuple{P,P}}) where {G, P<:Integer}
    # First create the standard RayTracingMesh
    standard_mesh = RayTracingMesh(faces, Ndiv)
    
    # Then convert to optimized version with spectral support
    optimMesh = RayTracingMeshOptim(standard_mesh)

    # Build spatial acceleration
    build_spatial_acceleration!(optimMesh)

    return optimMesh
end

function n_spectral_bins(mesh::RayTracingMeshOptim)
    return mesh.n_spectral_bins
end

function get_F_matrix(mesh::RayTracingMeshOptim, spectral_bin::Int=1)
    if mesh.spectral_mode == :spectral_variable
        return mesh.F_smooth[spectral_bin]
    else
        # For :grey and :spectral_uniform, use the single matrix
        return mesh.F_smooth
    end
end

function get_F_matrices(mesh::RayTracingMeshOptim)
    if mesh.spectral_mode == :spectral_variable
        return mesh.F_smooth
    else
        # For :grey and :spectral_uniform, return single matrix as vector for consistent interface
        return [mesh.F_smooth]
    end
end

# Quadrilateral constructor with spectral support - FIXED
function PolyFace2D{G}(p::SVector{4,Point2{G}}, b::SVector{4,Bool}, 
                         n_spectral_bins::Int=1,
                         kappa_default::G=zero(G), sigma_s_default::G=zero(G)) where {G}
    
    # Calculate geometric properties
    vertices = [p[1], p[2], p[3], p[4]]
    solidWalls = [b[1], b[2], b[3], b[4]]
    midPoint = convert(Point2{G}, (p[1]+p[2]+p[3]+p[4])/4)
    wallMidPoints = [convert(Point2{G}, (p[1]+p[2])/2), 
                     convert(Point2{G}, (p[2]+p[3])/2), 
                     convert(Point2{G}, (p[3]+p[4])/2), 
                     convert(Point2{G}, (p[4]+p[1])/2)]
    
    outwardNormals = [calculate_normal(p[1], p[2], midPoint),
                      calculate_normal(p[2], p[3], midPoint),
                      calculate_normal(p[3], p[4], midPoint),
                      calculate_normal(p[4], p[1], midPoint)]
    
    volume = convert(G, 0.5)*(p[1][1]*(p[2][2]-p[3][2])+p[2][1]*(p[3][2]-p[1][2])+p[3][1]*(p[1][2]-p[2][2])) + 
             convert(G, 0.5)*(p[3][1]*(p[4][2]-p[1][2])+p[4][1]*(p[1][2]-p[3][2])+p[1][1]*(p[3][2]-p[4][2]))
    
    area = [convert(G, norm(p[i]-p[mod(i, 4)+1])) for i = 1:4]
    subFaces = PolyFace2D{G}[]
    
    # Initialize extinction properties based on spectral mode
    if n_spectral_bins == 1
        kappa_g = kappa_default
        sigma_s_g = sigma_s_default
        epsilon = zeros(G, 4)
        
        # Grey mode - scalar values for volume properties
        j_g = zero(G)
        g_a_g = zero(G)
        e_g = zero(G)
        r_g = zero(G)
        g_g = zero(G)
        i_g = zero(G)
        q_in_g = zero(G)
        q_g = zero(G)
        T_in_g = zero(G)
        T_g = zero(G)
        
        # Grey mode - scalar values for wall properties  
        j_w = [zero(G), zero(G), zero(G), zero(G)]
        g_a_w = [zero(G), zero(G), zero(G), zero(G)]
        e_w = [zero(G), zero(G), zero(G), zero(G)]
        r_w = [zero(G), zero(G), zero(G), zero(G)]
        g_w = [zero(G), zero(G), zero(G), zero(G)]
        i_w = [zero(G), zero(G), zero(G), zero(G)]
        q_in_w = [zero(G), zero(G), zero(G), zero(G)]
        q_w = [zero(G), zero(G), zero(G), zero(G)]
        T_in_w = [zero(G), zero(G), zero(G), zero(G)]
        T_w = [zero(G), zero(G), zero(G), zero(G)]
    else
        kappa_g = fill(kappa_default, n_spectral_bins)
        sigma_s_g = fill(sigma_s_default, n_spectral_bins)
        epsilon = [zeros(G, n_spectral_bins) for _ in 1:4]
    
        # Spectral mode - vector values for volume properties
        j_g = zeros(G, n_spectral_bins)
        g_a_g = zeros(G, n_spectral_bins)
        e_g = zeros(G, n_spectral_bins)
        r_g = zeros(G, n_spectral_bins)
        g_g = zeros(G, n_spectral_bins)
        i_g = zeros(G, n_spectral_bins)
        q_in_g = zero(G)
        q_g = zero(G)
        T_in_g = zero(G)
        T_g = zero(G)
        
        # Spectral mode - vector values for wall properties
        j_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        g_a_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        e_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        r_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        g_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        i_w = [zeros(G, n_spectral_bins) for _ in 1:4]
        q_in_w = [zero(G) for _ in 1:4]
        q_w = [zero(G) for _ in 1:4]
        T_in_w = [zero(G) for _ in 1:4]
        T_w = [zero(G) for _ in 1:4]
    end
    
    # Call the default constructor with all fields
    return PolyFace2D{G}(
        vertices, solidWalls, midPoint, wallMidPoints, outwardNormals,
        volume, area, subFaces, epsilon,
        kappa_g, sigma_s_g,
        j_g, g_a_g, e_g, r_g, g_g, i_g, q_in_g, q_g, T_in_g, T_g,
        j_w, g_a_w, e_w, r_w, g_w, i_w, q_in_w, q_w, T_in_w, T_w
    )
end

# Triangle constructor with spectral support - FIXED
function PolyFace2D{G}(p::SVector{3, Point2{G}}, b::SVector{3,Bool}, 
                         n_spectral_bins::Int=1,
                         kappa_default::G=zero(G), sigma_s_default::G=zero(G)) where {G}
    
    # Calculate geometric properties
    vertices = [p[1], p[2], p[3]]
    solidWalls = [b[1], b[2], b[3]]
    midPoint = convert(Point2{G}, (p[1]+p[2]+p[3])./3)
    wallMidPoints = [convert(Point2{G}, (p[1]+p[2])/2), 
                     convert(Point2{G}, (p[2]+p[3])/2), 
                     convert(Point2{G}, (p[3]+p[1])/2)]
    
    outwardNormals = [calculate_normal(p[1], p[2], midPoint),
                      calculate_normal(p[2], p[3], midPoint),
                      calculate_normal(p[3], p[1], midPoint)]
    
    volume = convert(G, 0.5)*(p[1][1]*(p[2][2]-p[3][2])+p[2][1]*(p[3][2]-p[1][2])+p[3][1]*(p[1][2]-p[2][2]))
    
    area = [convert(G, norm(p[i]-p[mod(i, 3)+1])) for i = 1:3]
    subFaces = PolyFace2D{G}[]
    
    # Initialize extinction properties based on spectral mode
    if n_spectral_bins == 1
        kappa_g = kappa_default
        sigma_s_g = sigma_s_default
        epsilon = zeros(G, 3)
        
        # Grey mode - scalar values for volume properties
        j_g = zero(G)
        g_a_g = zero(G)
        e_g = zero(G)
        r_g = zero(G)
        g_g = zero(G)
        i_g = zero(G)
        q_in_g = zero(G)
        q_g = zero(G)
        T_in_g = zero(G)
        T_g = zero(G)
        
        # Grey mode - scalar values for wall properties  
        j_w = [zero(G), zero(G), zero(G)]
        g_a_w = [zero(G), zero(G), zero(G)]
        e_w = [zero(G), zero(G), zero(G)]
        r_w = [zero(G), zero(G), zero(G)]
        g_w = [zero(G), zero(G), zero(G)]
        i_w = [zero(G), zero(G), zero(G)]
        q_in_w = [zero(G), zero(G), zero(G)]
        q_w = [zero(G), zero(G), zero(G)]
        T_in_w = [zero(G), zero(G), zero(G)]
        T_w = [zero(G), zero(G), zero(G)]
    else
        kappa_g = fill(kappa_default, n_spectral_bins)
        sigma_s_g = fill(sigma_s_default, n_spectral_bins)
        epsilon = [zeros(G, n_spectral_bins) for _ in 1:3]
    
        # Spectral mode - vector values for volume properties
        j_g = zeros(G, n_spectral_bins)
        g_a_g = zeros(G, n_spectral_bins)
        e_g = zeros(G, n_spectral_bins)
        r_g = zeros(G, n_spectral_bins)
        g_g = zeros(G, n_spectral_bins)
        i_g = zeros(G, n_spectral_bins)
        q_in_g = zero(G)
        q_g = zero(G)
        T_in_g = zero(G)
        T_g = zero(G)
        
        # Spectral mode - vector values for wall properties
        j_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        g_a_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        e_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        r_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        g_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        i_w = [zeros(G, n_spectral_bins) for _ in 1:3]
        q_in_w = [zero(G) for _ in 1:3]
        q_w = [zero(G) for _ in 1:3]
        T_in_w = [zero(G) for _ in 1:3]
        T_w = [zero(G) for _ in 1:3]
    end
    
    # Call the default constructor with all fields
    return PolyFace2D{G}(
        vertices, solidWalls, midPoint, wallMidPoints, outwardNormals,
        volume, area, subFaces, epsilon,
        kappa_g, sigma_s_g,
        j_g, g_a_g, e_g, r_g, g_g, i_g, q_in_g, q_g, T_in_g, T_g,
        j_w, g_a_w, e_w, r_w, g_w, i_w, q_in_w, q_w, T_in_w, T_w
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
function build_uniform_grid(faces::Vector{PolyFace2D{G}}, cell_size::G) where {G}
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
function compute_bounding_boxes_opt(faces::Vector{PolyFace2D{G}}) where {G}
    return [BoundingBox(
        minimum(v[1] for v in face.vertices),
        maximum(v[1] for v in face.vertices),
        minimum(v[2] for v in face.vertices),
        maximum(v[2] for v in face.vertices)
    ) for face in faces]
end

# Build optimized spatial structure for a mesh
function build_optimized_spatial_structure(faces::Vector{PolyFace2D{G}}) where {G}
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

function build_spatial_grid(mesh::Vector{PolyFace2D{G}}, ncells::P) where {G,P<:Integer}
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
    Checks if extinction properties are uniform across all faces and spectral bins.
    Returns true if all faces have approximately the same kappa_g and sigma_s_g values
    for all spectral bins, false otherwise. Used to set global uniform_extinction flag.
    """
    
    first_beta = nothing
    first_n_bins = nothing
    is_spectral_mesh = false
    
    for (coarse_idx, coarse_face) in enumerate(mesh.coarse_mesh)
        for (fine_idx, fine_face) in enumerate(mesh.fine_mesh[coarse_idx])
            
            # Determine if this face is spectral
            current_is_spectral = isa(fine_face.kappa_g, Vector)
            
            # Check consistency of spectral structure across mesh
            if first_n_bins === nothing
                first_n_bins = current_is_spectral ? length(fine_face.kappa_g) : 1
                is_spectral_mesh = current_is_spectral
            else
                current_n_bins = current_is_spectral ? length(fine_face.kappa_g) : 1
                if current_n_bins != first_n_bins || current_is_spectral != is_spectral_mesh
                    error("Inconsistent spectral structure across mesh faces")
                end
            end
            
            # Check extinction values for uniformity
            if current_is_spectral
                # Spectral case - check each bin
                for bin in 1:length(fine_face.kappa_g)
                    current_beta = fine_face.kappa_g[bin] + fine_face.sigma_s_g[bin]
                    
                    if first_beta === nothing
                        first_beta = current_beta
                    else
                        if !isapprox(current_beta, first_beta, rtol=rtol)
                            println("Variable extinction detected - using variable ray tracing")
                            println("  Found spectral variation: bin $bin, face ($coarse_idx,$fine_idx)")
                            println("  β = $current_beta vs reference β = $first_beta")
                            return false
                        end
                    end
                end
            else
                # Grey case - single value
                current_beta = fine_face.kappa_g + fine_face.sigma_s_g
                
                if first_beta === nothing
                    first_beta = current_beta
                else
                    if !isapprox(current_beta, first_beta, rtol=rtol)
                        println("Variable extinction detected - using variable ray tracing")
                        println("  Found spatial variation: face ($coarse_idx,$fine_idx)")
                        println("  β = $current_beta vs reference β = $first_beta")
                        return false
                    end
                end
            end
        end
    end
    
    if mesh.spectral_mode != :grey
        println("Uniform extinction detected across $first_n_bins spectral bins (β_g=$(first_beta)) - using fast uniform ray tracing")
    else
        println("Uniform extinction detected (β_g=$(first_beta)) - using fast uniform ray tracing")
    end
    
    return true
end