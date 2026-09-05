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
    F_raw::Union{AbstractMatrix, SparseMatrixCSC{Float64, Int64}, Vector{SparseMatrixCSC{Float64, Int64}},
            Vector{AbstractMatrix}, Vector{AbstractMatrix{Float64}}, Vector{Matrix{Float64}}}              # UPDATED: Union for spectral support
    F_smooth::Union{AbstractMatrix, SparseMatrixCSC{Float64, Int64}, Vector{SparseMatrixCSC{Float64, Int64}},
            Vector{AbstractMatrix}, Vector{AbstractMatrix{Float64}}, Vector{Matrix{Float64}}}           # UPDATED: Union for spectral support
    surface_areas::VT
    volumes::VT
    surface_mapping::DIII
    volume_mapping::DII
    
    # NEW: Spectral metadata
    spectral_mode::Symbol        # :grey, :spectral_uniform, :spectral_variable
    n_spectral_bins::Int        # Number of spectral bins (1 for grey)
    wavelength_band_limits::Union{Nothing, Vector{Float64}}  # Wavelength boundaries [μm]
    surfaces_only::Bool         # indicates if the mesh includes volumes
    uniform_across_bin::Vector{Float64} # vector of uniform extinction (-1.0 where nonuniform)
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

mutable struct ViewFactorDomain3D{G,P<:Integer} <: SurfaceDomain3D{G,P}
    points::AbstractMatrix
    faces::AbstractMatrix
    Ndims::P
    facesMesh::Vector{PolyFace3D{G}}
    F_raw::AbstractMatrix
    F_smooth::AbstractMatrix  # View factors (wavelength-independent!)
    
    # NEW: Spectral metadata
    spectral_mode::Symbol        # :grey or :spectral
    n_spectral_bins::Int        # Number of spectral bins (1 for grey)
    wavelength_band_limits::Union{Nothing, Vector{G}}  # Wavelength boundaries [μm]
    energy_error::Union{Nothing, G, Vector{G}}
    uniform_epsilon::Bool        # Whether to use uniform epsilon solver
    surfaces_only::Bool         # dummy, always true, used for dispatch
end

# --- for recording rays ---
struct RayRecorder{P}
    ids::Vector{P}
    bin::P
    origins::Vector{Vector{Point2{Float64}}}
    endpoints::Vector{Vector{Point2{Float64}}}
end

# ---------- macro-surface grouping ----------
struct MacroSurface
    normal::Tuple{Float64,Float64}   # unit inward normal (orientation matters)
    offset::Float64                  # n . x for points on the line
    weight::Float64                  # summed element areas
    endpoints::Vector{NTuple{2,Float64}}   # element endpoints (for tests)
    element_ids::Vector{Int}         # global exchange-matrix indices
end

###################################################

"""
Intersection primitive. Every subface is triangulated for intersection;
a quad becomes two of these, both carrying the same `element`.
Edges precomputed for Möller–Trumbore.
"""
struct Tri3D{T<:AbstractFloat}
    a::Point3{T}
    ab::Vec3{T}
    ac::Vec3{T}
    element::Int32          # row/column index in F
end

"""
One radiative element, flattened. This is a row of F.
Triangles store v[4] == v[3].
"""
struct Facet3D{T<:AbstractFloat}
    v::SVector{4,Point3{T}}
    nv::Int8                # 3 or 4
    inwardNormal::Vec3{T}         # unit, into the enclosure
    e1::Vec3{T}             # tangent frame for cosine-weighted emission
    e2::Vec3{T}
    centroid::Point3{T}
    area::T
    tri_split::T            # area fraction of first triangle (quad sampling)
    parent::Int32           # index into facesMesh
end

struct BVHNode{T<:AbstractFloat}
    bmin::Point3{T}
    bmax::Point3{T}
    left::Int32             # left child; right is left+1. <0 marks a leaf
    first::Int32            # first primitive in `order` (leaves only)
    count::Int32
end

struct BVH{T<:AbstractFloat}
    nodes::Vector{BVHNode{T}}
    order::Vector{Int32}    # permutation into domain.tris
end

mutable struct RayTracingDomain3D_surfaces{G,P<:Integer} <: SurfaceDomain3D{G,P}
    points::AbstractMatrix
    faces::AbstractMatrix
    Ndims::P
    facesMesh::Vector{PolyFace3D{G}}

    F_raw::AbstractMatrix           # sparse, wavelength-independent
    F_smooth::AbstractMatrix

    # derived cache, rebuilt by flatten!
    facets::Union{Nothing, Vector{Facet3D{G}}}
    tris::Union{Nothing, Vector{Tri3D{G}}}
    bvh::Union{Nothing, BVH{G}}
    flattened::Bool

    spectral_mode::Symbol
    n_spectral_bins::Int
    wavelength_band_limits::Union{Nothing, Vector{G}}
    energy_error::Union{Nothing, G, Vector{G}}
    uniform_epsilon::Bool
    surfaces_only::Bool             # always true; kept for compatibility
end

struct DomainPlot3D{P,S,D,F}
    plot::P
    selected::S
    describe::D
    elements::F
    parents::Vector{Int}
end