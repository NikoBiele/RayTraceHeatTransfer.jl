#### functions to project a plane of 2d points back and forth between different planes

function project_plane_flat(plane_points)
    # Ensure we have at least 3 non-collinear points
    @assert length(plane_points) >= 3 "Need at least 3 points to define a plane"
    
    # Calculate the normal vector of the plane
    v1 = plane_points[2] - plane_points[1]
    v2 = plane_points[3] - plane_points[1]
    normal = normalize(cross(v1, v2))
    
    # Calculate rotation to align normal with z-axis
    z_axis = [0, 0, 1]
    rotation_axis = cross(normal, z_axis)
    rotation_angle = acos(clamp(dot(normal, z_axis), -1, 1))
    
    if norm(rotation_axis) < 1e-10
        # Normal is already aligned with z-axis (or its opposite)
        if dot(normal, z_axis) > 0
            R = Matrix{Float64}(I, 3, 3)
        else
            # Flip the plane
            R = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0]
        end
    else
        # Use quaternion rotation for better numerical stability
        rotation_axis = normalize(rotation_axis)
        q = [cos(rotation_angle/2); sin(rotation_angle/2) * rotation_axis...]
        R = quat_to_rot_matrix(q)
    end
    
    # Calculate translation to move a point on the plane to the origin
    T = -plane_points[1]
    
    return R, T
end

function quat_to_rot_matrix(q)
    w, x, y, z = q
    return [
        1-2y^2-2z^2  2x*y-2z*w    2x*z+2y*w;
        2x*y+2z*w    1-2x^2-2z^2  2y*z-2x*w;
        2x*z-2y*w    2y*z+2x*w    1-2x^2-2y^2
    ]
end

function project_point_flat(point, R, T)
    return R * (point + T)
end

function project_point_back(xy_point, R, T)
    return R \ xy_point - T
end

function project_face_flat(cube, face)
    face_points = [cube[i, :] for i in face]
    R, T = project_plane_flat(face_points)
    projected_points = [project_point_flat(point, R, T) for point in face_points]
    return projected_points, R, T
end

function project_face_back(flat_cube_points, R, T)
    num_faces = length(flat_cube_points)
    vector_length = size(flat_cube_points[1])
    projected_points = Vector{Vector{Vector{Float64}}}(undef, 4)
    for i = 1:4
        projected_points[i] = Vector{Vector{Float64}}(undef, vector_length[1])
        for j = 1:vector_length[1]
            projected_points[i][j] = project_point_back(flat_cube_points[i][j], R, T)
        end
    end
    return projected_points
end