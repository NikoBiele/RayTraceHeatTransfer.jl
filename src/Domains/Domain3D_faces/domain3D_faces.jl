mutable struct Face3D{G}

    # geometric variables (unchanged)
    vertices::Vector{Point3{G}}
    midPoint::Point3{G}
    inwardNormal::Point3{G}
    area::G
    subFaces::Union{Nothing, Vector{Face3D{G}}}

    # UPDATED: emissivity - now Union for spectral support
    epsilon::Union{G, Vector{G}}  # scalar (grey) or vector (spectral)

    # UPDATED: state variables - spectral quantities are vectors, physical quantities are scalars
    j_w::Union{Nothing, G, Vector{G}}     # outgoing power [W] - spectral
    g_a_w::Union{Nothing, G, Vector{G}}   # incident absorbed power [W] - spectral
    e_w::Union{Nothing, G, Vector{G}}     # emissive power [W] - spectral
    r_w::Union{Nothing, G, Vector{G}}     # reflected power [W] - spectral
    g_w::Union{Nothing, G, Vector{G}}     # incident power [W] - spectral
    i_w::Union{Nothing, G, Vector{G}}     # total intensity [W*m^(-2)*sr^(-1)] - spectral
    q_in_w::G          # input source term [W] - scalar (total)
    q_w::Union{Nothing, G}                # source term [W] - scalar (total)
    T_in_w::G          # input temperature [K] - scalar
    T_w::Union{Nothing, G}                # temperature [K] - scalar (physical quantity)
end

# Updated constructor with spectral support
function Face3D(vertices::Vector{Point3{G}}, domain_midpoint::Point3{G}, 
                epsilon::Union{G, Vector{G}}, q_in_w::G, T_in_w::G;
                n_spectral_bins::Int=1) where {G}
    
    midPoint = sum(vertices)/length(vertices)
    
    if length(vertices) == 3
        area = norm(cross(vertices[2] - vertices[1], vertices[3] - vertices[1]))/2
    else
        area = norm(cross(vertices[2] - vertices[1], vertices[4] - vertices[1]))
    end
    
    inwardNormal = calculate_normal(vertices[1], vertices[2], vertices[3], domain_midpoint)
    
    # Initialize spectral properties based on epsilon type
    is_spectral = isa(epsilon, Vector)
    
    if is_spectral
        n_bins = length(epsilon)
        # Spectral mode - radiative quantities are vectors
        j_w = zeros(G, n_bins)
        g_a_w = zeros(G, n_bins)
        e_w = zeros(G, n_bins)
        r_w = zeros(G, n_bins)
        g_w = zeros(G, n_bins)
        i_w = zeros(G, n_bins)
        # Physical quantities remain scalar
        q_w = nothing
        T_w = nothing
    else
        # Grey mode - all scalar values
        j_w = nothing
        g_a_w = nothing
        e_w = nothing
        r_w = nothing
        g_w = nothing
        i_w = nothing
        q_w = nothing
        T_w = nothing
    end
    
    return Face3D(vertices, midPoint, inwardNormal, area, nothing, epsilon, 
                  j_w, g_a_w, e_w, r_w, g_w, i_w,
                  q_in_w, q_w, T_in_w, T_w)
end

# Helper function unchanged
function calculate_normal(p1::Point3{G}, p2::Point3{G}, p3::Point3{G}, midpoint::Point3{G}) where {G}
    edge1 = p2 - p1
    edge2 = p3 - p1
    normal = normalize(cross(edge1, edge2))
    
    face_midpoint = (p1 + p2 + p3) / 3
    if dot(normal, midpoint - face_midpoint) < 0
        normal = -normal
    end
    
    return normal
end

mutable struct Domain3D_faces{G,P<:Integer}
    points::Matrix{G}
    faces::Matrix{P}
    Ndims::P
    facesMesh::Vector{Face3D{G}}
    F::Matrix{G}  # View factors (wavelength-independent!)
    
    # NEW: Spectral metadata
    spectral_mode::Symbol        # :grey or :spectral
    n_spectral_bins::Int        # Number of spectral bins (1 for grey)
    wavelength_band_limits::Union{Nothing, Vector{G}}  # Wavelength boundaries [μm]
    energy_error::Union{Nothing, G, Vector{G}}
end

# Updated constructor with spectral support
function Domain3D_faces(points::Matrix{G}, faces::Matrix{P}, Ndims::P,
                       q_in_w::Vector{G}, T_in_w::Vector{G}, 
                       epsilon::Union{Vector{G}, Vector{Vector{G}}}) where {G, P<:Integer}
    
    # Determine spectral mode from epsilon
    is_spectral = isa(epsilon[1], Vector)
    if is_spectral
        if std([std(epsilon[i]) for i in 1:size(epsilon, 1)]) > 1e-6
            spectral_mode = :spectral_variable
        else
            spectral_mode = :spectral_uniform
        end
    else
        spectral_mode = :grey
    end
    n_bins = is_spectral ? length(epsilon[1]) : 1
    
    # Calculate domain midpoint
    domain_mid = sum(points, dims=1)/size(points, 1)    
    domain_midpoint = Point3{G}(domain_mid[1], domain_mid[2], domain_mid[3])
    
    # Create super faces
    superFaces = Face3D[]
    for (i, face_rows) in enumerate(eachrow(faces))
        if length(face_rows) == 4
            points3d = [Point3{G}(points[face_rows[j],:]) for j in 1:4]
        else
            points3d = [Point3{G}(points[face_rows[j],:]) for j in 1:3]
        end
        
        push!(superFaces, Face3D(points3d, domain_midpoint, epsilon[i], q_in_w[i], T_in_w[i]))
    end
    
    # Mesh the faces
    mesh3D = mesh3D_faces(points, faces, Ndims)
    num_faces = length(mesh3D)
    num_points = length(mesh3D[1][1])
    
    for i in 1:num_faces
        superFaces[i].subFaces = Face3D[]
        
        # First pass: create all subfaces with q=0
        for j in 1:num_points
            p1 = Point3{G}(mesh3D[i][1][j]...)
            p2 = Point3{G}(mesh3D[i][2][j]...)
            p3 = Point3{G}(mesh3D[i][3][j]...)
            p4 = Point3{G}(mesh3D[i][4][j]...)
            
            if isapprox(p3, p4, atol=1e-5) 
                push!(superFaces[i].subFaces, 
                    Face3D([p1, p2, p3], domain_midpoint, epsilon[i], 0.0, T_in_w[i]))
            else
                push!(superFaces[i].subFaces, 
                    Face3D([p1, p2, p3, p4], domain_midpoint, epsilon[i], 0.0, T_in_w[i]))
            end
        end
        
        # Second pass: distribute flux proportional to area
        total_area = sum(subface.area for subface in superFaces[i].subFaces)
        
        for subface in superFaces[i].subFaces
            # Each subface gets flux proportional to its area fraction
            subface.q_in_w = q_in_w[i] * (subface.area / total_area)
        end
    end
    
    # Calculate view factors (wavelength-independent!)
    println("Computing view factors (geometry only, wavelength-independent)...")
    F = exchangeFactors3D!(superFaces)
    
    # Return with superFaces as facesMesh (this is the variable name expected)
    return Domain3D_faces{G, P}(points, faces, Ndims, superFaces, F, 
                                spectral_mode, n_bins, nothing, nothing)
end