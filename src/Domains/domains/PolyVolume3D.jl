# Primary constructor: Extrude a 2D PolyVolume2D into 3D with specified depth
# The 2D polygon lies in the x-y plane at z=0 (front face) and is extruded in the +z direction (into the page)
# Face ordering: side faces (one per 2D edge), then front face (z=0), then back face (z=depth)
# Front and back faces are ALWAYS solid
# Side faces ALWAYS inherit properties from the 2D walls
function PolyVolume3D{G}(poly2d::PolyVolume2D{G}, 
                         depth::G,
                         epsilon_front::Union{G, Vector{G}},
                         epsilon_back::Union{G, Vector{G}},
                         q_in_front::G,
                         q_in_back::G,
                         T_in_front::G,
                         T_in_back::G) where {G}
    
    n_edges = length(poly2d.vertices)
    
    # Side faces inherit from 2D walls
    solidSides = copy(poly2d.solidWalls)
    epsilon_sides = copy(poly2d.epsilon)
    q_in_sides = copy(poly2d.q_in_w)
    T_in_sides = copy(poly2d.T_in_w)
    
    # Determine spectral mode
    is_spectral = isa(epsilon_front, Vector)
    n_spectral_bins = is_spectral ? length(epsilon_front) : 1
    
    # Extinction coefficients
    kappa_g = poly2d.kappa_g
    sigma_s_g = poly2d.sigma_s_g
    
    # Create 3D vertices by extruding 2D vertices
    # Front face at z_start, back face at z_start+depth
    vertices_front = [Point3{G}(v[1], v[2], 0.0) for v in poly2d.vertices]
    vertices_back = [Point3{G}(v[1], v[2], depth) for v in poly2d.vertices]
    
    # All 3D vertices: front face vertices (1:n), then back face vertices (n+1:2n)
    vertices = vcat(vertices_front, vertices_back)
    
    # Calculate volume centroid
    midPoint = convert(Point3{G}, sum(vertices) / length(vertices))
    
    # Create face structures
    # Face ordering: side faces (1:n_edges), front face (n_edges+1), back face (n_edges+2)
    faces = Vector{PolyFace3D{G}}(undef, n_edges + 2)
    
    # Create side faces (rectangular faces formed by extruding each 2D edge)
    for i in 1:n_edges
        j = mod1(i + 1, n_edges)  # Next vertex index
        
        # Vertices for this rectangular side face (counter-clockwise from outside)
        # front_i -> front_j -> back_j -> back_i
        face_verts = [
            vertices_front[i],
            vertices_front[j],
            vertices_back[j],
            vertices_back[i]
        ]

        faces[i] = PolyFace3D(face_verts, solidSides[i], midPoint, epsilon_sides[i], q_in_sides[i], T_in_sides[i],
                              n_spectral_bins=n_spectral_bins)
    end

    # Create front face (z=0) - always solid
    front_verts = vertices_front
    faces[n_edges + 1] = PolyFace3D(front_verts, true, midPoint, epsilon_front, q_in_front, T_in_front,
                                     n_spectral_bins=n_spectral_bins)

    # Create back face (z=depth) - always solid
    back_verts = vertices_back
    faces[n_edges + 2] = PolyFace3D(back_verts, true, midPoint, epsilon_back, q_in_back, T_in_back,
                                     n_spectral_bins=n_spectral_bins)
 
    # Calculate volume: area of 2D polygon × depth
    volume = poly2d.volume * depth
    
    # Initialize empty subvolume list
    subVolumes = PolyVolume3D{G}[]
    
    # Initialize extinction properties based on spectral mode
    if n_spectral_bins == 1
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
    else        
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
    end
    
    # Call the default constructor with all fields
    return PolyVolume3D{G}(
        vertices, faces, midPoint, volume, subVolumes,
        kappa_g, sigma_s_g,
        j_g, g_a_g, e_g, r_g, g_g, i_g, q_in_g, q_g, T_in_g, T_g
    )
end

# Hexahedron (box) constructor with spectral support
# Vertices ordered: back face (1-4), front face (5-8)
# This is the low-level constructor for when you have explicit vertices
#   5---8
#  /|  /|
# 6---7 |
# | 1-|-4
# |/  |/
# 2---3
function PolyVolume3D{G}(vertices::Vector{Point3{G}}, 
                         solidFaces::Vector{Bool},
                         epsilon::Union{Vector{G}, Vector{Vector{G}}},  # emissivities for all faces
                         q_in_faces::Vector{G},  # input heat sources for all faces
                         T_in_faces::Vector{G},  # input temperatures for all faces
                         n_spectral_bins::Int=1,
                         kappa_default::G=zero(G), 
                         sigma_s_default::G=zero(G)) where {G}
    
    @assert length(vertices) == 8 "Hexahedron requires 8 vertices"
    n_faces = length(solidFaces)
    @assert length(solidFaces) == n_faces "Hexahedron has 6 faces"
    @assert length(epsilon) == n_faces "Must provide emissivity for all faces"
    @assert length(q_in_faces) == n_faces "Must provide q_in for all faces"
    @assert length(T_in_faces) == n_faces "Must provide T_in for all faces"
    
    # Calculate volume centroid
    midPoint = convert(Point3{G}, sum(vertices) / 8)
    
    # Define faces by vertex indices (counter-clockwise when viewed from outside)
    # Order: 4 side faces, back face, front face
    face_indices = [
        [1, 2, 6, 5],  # Side 1
        [2, 3, 7, 6],  # Side 2
        [3, 4, 8, 7],  # Side 3
        [4, 1, 5, 8],  # Side 4
        [4, 3, 2, 1],  # Back face (z=0, reversed)
        [5, 6, 7, 8]   # Front face (z=depth)
    ]
    
    # Create PolyFace3D objects for each face
    faces = Vector{PolyFace3D{G}}(undef, n_faces)
    for i in 1:n_faces
        face_verts = [vertices[j] for j in face_indices[i]]
        faces[i] = PolyFace3D(face_verts, solidFaces[i], midPoint, epsilon[i], q_in_faces[i], T_in_faces[i],
                              n_spectral_bins=n_spectral_bins)
    end
    
    # Calculate volume using divergence theorem or decomposition
    # For a box: V = ||(v7 - v1) · ((v3 - v1) × (v5 - v1))||
    volume = abs(dot(vertices[7] - vertices[1], 
                     cross(vertices[3] - vertices[1], vertices[5] - vertices[1])))
    
    # Initialize empty subvolume list
    subVolumes = PolyVolume3D{G}[]
    
    # Initialize extinction properties based on spectral mode
    if n_spectral_bins == 1
        # Grey mode
        kappa_g = kappa_default
        sigma_s_g = sigma_s_default
        
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
    else
        # Spectral mode
        kappa_g = fill(kappa_default, n_spectral_bins)
        sigma_s_g = fill(sigma_s_default, n_spectral_bins)
        
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
    end
    
    # Call the default constructor with all fields
    return PolyVolume3D{G}(
        vertices, faces, midPoint, volume, subVolumes,
        kappa_g, sigma_s_g,
        j_g, g_a_g, e_g, r_g, g_g, i_g, q_in_g, q_g, T_in_g, T_g
    )
end