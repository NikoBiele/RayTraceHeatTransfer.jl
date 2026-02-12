# Define a custom PolyVolume2D type
mutable struct PolyVolume2D{G}

    # geometric variables, fixed types, no uncertainty (unchanged)
    vertices::Vector{Point2{G}}
    solidWalls::Vector{Bool}
    midPoint::Point2{G}
    wallMidPoints::Vector{Point2{G}}
    inwardNormals::Vector{Point2{G}}
    volume::G
    area::Vector{G}
    subVolumes::Vector{PolyVolume2D{G}} # not fixed

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

struct GridCell2D{G,P}
    face_indices::Vector{P}
    bounds::Tuple{Point2{G}, Point2{G}}  # min and max points
end

struct SpatialGrid{G,P}
    cells::Matrix{GridCell2D{G,P}}
    cell_size::Point2{G}
    min_point::Point2{G}
    max_point::Point2{G}
    dims::Tuple{P,P}
end

mutable struct IntermediateMesh2D{VPF,VVPF,MT,GRID}
    coarse_mesh::VPF
    fine_mesh::VVPF
    coarse_grid::GRID  
    fine_grids::Vector{GRID}
    
    # UPDATED: Exchange factor matrices - Union for spectral support
    F_raw::Union{MT, Vector{MT}}              # Single matrix (grey) or vector of matrices (spectral)
    F_smooth::Union{MT, Vector{MT}}           # Single matrix (grey) or vector of matrices (spectral)
end

# Bounding box for fast rejection tests
struct BoundingBox2D{G}
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

# Enhanced RayTracingDomain2D with spectral support
mutable struct RayTracingDomain2D{VPF,VVPF,MT,VT,DIII,DII,GRID}
    # Original fields (matching your exact RayTracingMesh structure)
    coarse_mesh::VPF
    fine_mesh::VVPF
    coarse_grid::GRID
    fine_grids::Vector{GRID}
    F_raw::Union{MT, Vector{MT}}              # UPDATED: Union for spectral support
    F_smooth::Union{MT, Vector{MT}}           # UPDATED: Union for spectral support
    surface_areas::VT
    volumes::VT
    surface_mapping::DIII
    volume_mapping::DII
    
    # NEW: Spectral metadata
    spectral_mode::Symbol        # :grey, :spectral_uniform, :spectral_variable
    n_spectral_bins::Int        # Number of spectral bins (1 for grey)
    wavelength_band_limits::Union{Nothing, Vector{Float64}}  # Wavelength boundaries [μm]
    uniform_extinction::Bool     # Whether to use uniform exchange factors across all wavelengths
    
    # Optimized cache structures (existing)
    coarse_face_cache::Vector{PolyVolume2D}  # Flattened for direct indexing
    fine_face_cache::Vector{Vector{PolyVolume2D}}  # Pre-allocated fine faces
    
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
    coarse_bboxes_opt::Union{Nothing, Vector{BoundingBox2D}}
    fine_grids_opt::Union{Nothing, Vector{UniformGrid}}
    fine_bboxes_opt::Union{Nothing, Vector{Vector{BoundingBox2D}}}

    # automatic calculation of global energy conservation error
    energy_error::Union{Nothing, G, Vector{G}} where {G}
end

mutable struct PolyFace3D{G}

    # geometric variables (unchanged)
    vertices::Union{Nothing, Vector{Point3{G}}}
    solidFace::Union{Nothing, Bool}
    midPoint::Union{Nothing, Point3{G}}
    inwardNormal::Union{Nothing, Point3{G}}
    area::Union{Nothing, G}
    subFaces::Union{Nothing, Vector{PolyFace3D{G}}}

    # UPDATED: emissivity - now Union for spectral support
    epsilon::Union{Nothing, G, Vector{G}}  # scalar (grey) or vector (spectral)

    # UPDATED: state variables - spectral quantities are vectors, physical quantities are scalars
    j_w::Union{Nothing, G, Vector{G}}     # outgoing power [W] - spectral
    g_a_w::Union{Nothing, G, Vector{G}}   # incident absorbed power [W] - spectral
    e_w::Union{Nothing, G, Vector{G}}     # emissive power [W] - spectral
    r_w::Union{Nothing, G, Vector{G}}     # reflected power [W] - spectral
    g_w::Union{Nothing, G, Vector{G}}     # incident power [W] - spectral
    i_w::Union{Nothing, G, Vector{G}}     # total intensity [W*m^(-2)*sr^(-1)] - spectral
    q_in_w::Union{Nothing, G}          # input source term [W] - scalar (total)
    q_w::Union{Nothing, G}                # source term [W] - scalar (total)
    T_in_w::Union{Nothing, G}          # input temperature [K] - scalar
    T_w::Union{Nothing, G}                # temperature [K] - scalar (physical quantity)
end

mutable struct ViewFactorDomain3D{G,P<:Integer}
    points::Matrix{G}
    faces::Matrix{P}
    Ndims::P
    facesMesh::Vector{PolyFace3D{G}}
    F_raw::Matrix{G}
    F_smooth::Matrix{G}  # View factors (wavelength-independent!)
    
    # NEW: Spectral metadata
    spectral_mode::Symbol        # :grey or :spectral
    n_spectral_bins::Int        # Number of spectral bins (1 for grey)
    wavelength_band_limits::Union{Nothing, Vector{G}}  # Wavelength boundaries [μm]
    energy_error::Union{Nothing, G, Vector{G}}
    uniform_epsilon::Bool        # Whether to use uniform epsilon solver
end

#######################

# PolyVolume3D - 3D volumetric element for participating media ray tracing
# This is the 3D equivalent of PolyVolume2D, representing a hexahedral (box) element
# Uses nested PolyFace3D structures for faces (consistent with existing architecture)

mutable struct PolyVolume3D{G}
    
    # Geometric variables - fixed types, no uncertainty
    vertices::Union{Nothing, Vector{Point3{G}}}          # 8 vertices for hexahedron (box)
    faces::Union{Nothing, Vector{PolyFace3D{G}}} # up to 6 faces (each is a PolyFace3D with solidFace field)
    midPoint::Union{Nothing, Point3{G}}                  # Volume centroid
    volume::Union{Nothing, G}                            # 3D volume [m^3]
    subVolumes::Union{Nothing, Vector{Vector{PolyVolume3D{G}}}} # Hierarchical subdivision (one depth vector per 2d volume)
    
    # Local gas extinction properties - scalar (grey) or vector (spectral)
    kappa_g::Union{Nothing, G, Vector{G}}         # absorption coefficient [m^-1]
    sigma_s_g::Union{Nothing, G, Vector{G}}       # scattering coefficient [m^-1]
    
    # State variables (volume) - radiative quantities
    j_g::Union{Nothing, G, Vector{G}}             # outgoing power [W]
    g_a_g::Union{Nothing, G, Vector{G}}           # incident absorbed power [W]
    e_g::Union{Nothing, G, Vector{G}}             # emissive power [W]
    r_g::Union{Nothing, G, Vector{G}}             # reflected power [W]
    g_g::Union{Nothing, G, Vector{G}}             # incident power [W]
    i_g::Union{Nothing, G, Vector{G}}             # total intensity [W*m^(-2)*sr^(-1)]
    q_in_g::Union{Nothing, G}                     # input source term [W]
    q_g::Union{Nothing, G}                        # source term [W]
    T_in_g::Union{Nothing, G}                     # input temperature [K]
    T_g::Union{Nothing, G}                        # temperature [K]
end

# RayTracingDomain3D - Main domain structure for 3D participating media ray tracing
# This is the 3D equivalent of RayTracingDomain2D

mutable struct RayTracingDomain3D{G}
    # Mesh hierarchy
    coarse_mesh::Vector{PolyVolume3D{G}}                   # Vector{PolyVolume3D{G}}
    # fine_mesh::Vector{Vector{Vector{PolyVolume3D{G}}}}    # Vector{Vector{Vector{PolyVolume3D{G}}}}
    
    # Exchange factors
    F_raw::Union{Matrix{G}, Vector{Matrix{G}}}       # Raw exchange factors
    F_smooth::Union{Matrix{G}, Vector{Matrix{G}}}    # Smoothed exchange factors
    
    # Geometric data
    surface_areas::Vector{G}           # Areas of all solid faces [m^2]
    # volumes::Vector{G}                 # Volumes of all cells [m^3]
    
    # Index mappings
    surface_mapping::Dict{Tuple{Int,Int,Int}, Int}  # (coarse, fine, depth, face) -> global surface index
    # volume_mapping::Dict{Tuple{Int,Int,Int}, Int}   # (coarse, fine) -> global volume index, not necessary
    
    # Spectral metadata
    spectral_mode::Symbol              # :grey, :spectral_uniform, :spectral_variable
    n_spectral_bins::Int               # Number of spectral bins (1 for grey)
    wavelength_band_limits::Union{Nothing, Vector{Float64}}  # Wavelength boundaries [μm]
    uniform_extinction::Bool           # Whether extinction is uniform across domain

    # use 2d accelerations since it is extruded 2d (search 2d, then 1d)

    # Optimized cache structures (existing)
    from2d_coarse_face_cache::Union{Nothing, Vector{PolyVolume2D{G}}} # Flattened for direct indexing
    from2d_fine_face_cache::Union{Nothing, Vector{Vector{PolyVolume2D{G}}}}  # Pre-allocated fine faces
    
    # Pre-computed geometric data for faster access
    from2d_coarse_wall_normals::Union{Nothing, Vector{Vector{Point2}}}  # Outward normals per face
    from2d_coarse_wall_midpoints::Union{Nothing, Vector{Vector{Point2}}}  # Wall midpoints per face
    from2d_fine_wall_normals::Union{Nothing, Vector{Vector{Vector{Point2}}}}  # Fine mesh normals
    from2d_fine_wall_midpoints::Union{Nothing, Vector{Vector{Vector{Point2}}}}  # Fine mesh wall midpoints
    
    # Spatial acceleration structures (existing)
    from2d_coarse_bounding_boxes::Union{Nothing, Vector{Tuple{Point2, Point2}}}  # (min, max) per coarse face
    from2d_fine_bounding_boxes::Union{Nothing, Vector{Vector{Tuple{Point2, Point2}}}}  # Bounding boxes for fine faces
    
    # Optimized spatial acceleration structures
    from2d_coarse_grid_opt::Union{Nothing, UniformGrid}
    from2d_coarse_bboxes_opt::Union{Nothing, Vector{BoundingBox2D}}
    from2d_fine_grids_opt::Union{Nothing, Vector{UniformGrid}}
    from2d_fine_bboxes_opt::Union{Nothing, Vector{Vector{BoundingBox2D}}}
    
    # Energy conservation tracking
    energy_error::Union{Nothing, G, Vector{G}}

end