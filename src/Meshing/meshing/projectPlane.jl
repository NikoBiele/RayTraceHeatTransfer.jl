#### functions to project a plane of 2d points back and forth between different planes
#### (all arithmetic in the element type of the input points)

function projectPlaneFlat(plane_points)
    # Ensure we have at least 3 non-collinear points
    # @assert length(plane_points) >= 3 "Need at least 3 points to define a plane"
    
    # Calculate the normal vector of the plane
    v1 = plane_points[2] - plane_points[1]
    v2 = plane_points[3] - plane_points[1]
    normal = normalize(cross(v1, v2))
    G = eltype(normal)                                   # working precision
    
    # Calculate rotation to align normal with z-axis
    z_axis = G[0, 0, 1]
    rotation_axis = cross(normal, z_axis)
    # angle between normal and z-axis; atan(|n × z|, n · z) is accurate for all
    # angles, whereas acos loses precision when the normal is nearly aligned
    rotation_angle = atan(norm(rotation_axis), dot(normal, z_axis))
    
    if norm(rotation_axis) < 100 * eps(G)
        # Normal is already aligned with z-axis (or its opposite), to working precision
        if dot(normal, z_axis) > 0
            R = Matrix{G}(I, 3, 3)
        else
            # Flip the plane
            R = G[1 0 0; 0 1 0; 0 0 -1]
        end
    else
        # Use quaternion rotation for better numerical stability
        rotation_axis = normalize(rotation_axis)
        q = [cos(rotation_angle/2); sin(rotation_angle/2) * rotation_axis...]
        R = quatToRotMatrix(q)
    end
    
    # Calculate translation to move a point on the plane to the origin
    T = -plane_points[1]
    
    return R, T
end

function quatToRotMatrix(q)
    w, x, y, z = q
    return [
        1-2y^2-2z^2  2x*y-2z*w    2x*z+2y*w;
        2x*y+2z*w    1-2x^2-2z^2  2y*z-2x*w;
        2x*z-2y*w    2y*z+2x*w    1-2x^2-2y^2
    ]
end

function projectPointFlat(point, R, T)
    return R * (point + T)
end

function projectPointBack(xy_point, R, T)
    return R \ xy_point - T
end

function projectFaceFlat(cube, face)
    face_points = [cube[i, :] for i in face]
    R, T = projectPlaneFlat(face_points)
    projected_points = [projectPointFlat(point, R, T) for point in face_points]
    return projected_points, R, T
end

function projectFaceBack(flat_cube_points, R, T)
    num_faces = length(flat_cube_points)
    vector_length = size(flat_cube_points[1])
    S = promote_type(eltype(R), eltype(T))               # keep the working precision
    projected_points = Vector{Vector{Vector{S}}}(undef, 4)
    for i = 1:4
        projected_points[i] = Vector{Vector{S}}(undef, vector_length[1])
        for j = 1:vector_length[1]
            projected_points[i][j] = projectPointBack(flat_cube_points[i][j], R, T)
        end
    end
    return projected_points
end