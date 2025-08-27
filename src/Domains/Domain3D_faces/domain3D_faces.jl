mutable struct Face3D{G,T}

    # geometric variables
    vertices::Vector{Point3{G}}
    midPoint::Point3{G}
    inwardNormal::Point3{G}
    area::G
    subFaces::Union{Nothing, Vector{Face3D{G,T}}}

    # variables with possibility of uncertainty

    # emissivity
    epsilon::T

    # state variables (walls)
    j_w::Union{Nothing, T} # outgoing power [W]
    g_a_w::Union{Nothing, T} # incident absorbed power [W]
    e_w::Union{Nothing, T} # emissive power [W]
    r_w::Union{Nothing, T} # reflected power [W]
    g_w::Union{Nothing, T} # incident power [W]
    i_w::Union{Nothing, T} # total intensity [W*m^(-2)*sr^(-1)]
    q_in_w::T # input source term [W]
    q_w::Union{Nothing, T} # source term [W]
    T_in_w::T # input temperature [K]
    T_w::Union{Nothing, T} # temperature [K]
end

function Face3D(vertices::Vector{Point3{G}}, domain_midpoint::Point3{G}, epsilon::G,
                q_in_w::G, T_in_w::G) where {G}
    midPoint = sum(vertices)/length(vertices)
    if length(vertices) == 3
        area = norm(cross(vertices[2] - vertices[1], vertices[3] - vertices[1]))/2
    else
        area = norm(cross(vertices[2] - vertices[1], vertices[4] - vertices[1]))
    end
    inwardNormal = calculate_normal(vertices[1], vertices[2], vertices[3], domain_midpoint)
    return Face3D(vertices, midPoint, inwardNormal, area, nothing, epsilon, 
                nothing, nothing, nothing, nothing, nothing, nothing,
                q_in_w, nothing, T_in_w, nothing)
end

function calculate_normal(p1::Point3{G}, p2::Point3{G}, p3::Point3{G}, midpoint::Point3{G}) where G
    # Calculate two edges of the triangle
    edge1 = p2 - p1
    edge2 = p3 - p1
    
    # Cross product gives the normal (right-hand rule)
    normal = normalize(cross(edge1, edge2))
    
    # Check if normal is pointing toward the domain midpoint
    face_midpoint = (p1 + p2 + p3) / 3
    if dot(normal, midpoint - face_midpoint) < 0
        normal = -normal  # Flip if pointing away from domain
    end
    
    return normal
end

mutable struct Domain3D_faces{P<:AbstractFloat,K<:Integer}
    points::Matrix{P}
    faces::Matrix{K}
    Ndims::K
    facesMesh::Vector{Face3D{P,P}}
    F::Matrix{P}
end

function Domain3D_faces(points::Matrix{P}, faces::Matrix{K}, Ndims::K,
                q_in_w::Vector{P}, T_in_w::Vector{P}, epsilon::Vector{P}) where {P<:AbstractFloat, K<:Integer}
    domain_mid = sum(points, dims=1)/size(points, 1)    
    domain_midpoint = Point3{P}(domain_mid[1], domain_mid[2], domain_mid[3])
    superFaces = Face3D[]
    for (i, face_rows) in enumerate(eachrow(faces))
        if length(face_rows) == 4
            points3d = [Point3{P}(points[face_rows[1],:]), Point3{P}(points[face_rows[2],:]), Point3{P}(points[face_rows[3],:]), Point3{P}(points[face_rows[4],:])]
        else
            points3d = [Point3{P}(points[face_rows[1],:]), Point3{P}(points[face_rows[2],:]), Point3{P}(points[face_rows[3],:])]
        end
        push!(superFaces, Face3D(points3d, domain_midpoint, epsilon[i], q_in_w[i], T_in_w[i]))
    end
    mesh3D = mesh3D_faces(points, faces, Ndims)
    num_faces = length(mesh3D)
    num_points = length(mesh3D[1][1])
    for i in 1:num_faces
        superFaces[i].subFaces = Face3D[]
        for j in 1:num_points
            p1 = Point3{P}(mesh3D[i][1][j][1], mesh3D[i][1][j][2], mesh3D[i][1][j][3])
            p2 = Point3{P}(mesh3D[i][2][j][1], mesh3D[i][2][j][2], mesh3D[i][2][j][3])
            p3 = Point3{P}(mesh3D[i][3][j][1], mesh3D[i][3][j][2], mesh3D[i][3][j][3])
            p4 = Point3{P}(mesh3D[i][4][j][1], mesh3D[i][4][j][2], mesh3D[i][4][j][3])
            if isapprox(p3, p4, atol=1e-5) 
                push!(superFaces[i].subFaces, Face3D([p1, p2, p3], domain_midpoint, epsilon[i], q_in_w[i], T_in_w[i]))
            else
                push!(superFaces[i].subFaces, Face3D([p1, p2, p3, p4], domain_midpoint, epsilon[i], q_in_w[i], T_in_w[i]))
            end
        end
    end            
    F = exchangeFactors3D!(superFaces)

    return Domain3D_faces{P, K}(points, faces, Ndims, superFaces, F)
end
