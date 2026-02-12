# Updated 2d version of meshTriangle with spectral support
function meshTriangle(face::PolyVolume2D{G}, Ndim::P) where {G, P<:Integer}
    # Determine spectral mode from parent face
    is_spectral = isa(face.kappa_g, Vector)
    n_bins = is_spectral ? length(face.kappa_g) : 1
    
    # Get default extinction values from parent face
    kappa_default = is_spectral ? face.kappa_g[1] : face.kappa_g
    sigma_s_default = is_spectral ? face.sigma_s_g[1] : face.sigma_s_g
    
    # save triangle midpoint
    triangle_midPoint = face.midPoint
    
    # insert a point opposite of one the point opposite the longest edge
    norm1 = norm(face.vertices[1] - face.vertices[2])
    norm2 = norm(face.vertices[2] - face.vertices[3])
    norm3 = norm(face.vertices[3] - face.vertices[1])
    max_norm, max_index = findmax([norm1, norm2, norm3])
    
    # mirror the opposite point
    if max_index == 1
        point_to_mirror = face.vertices[3]
        startpoint = face.vertices[1]
        line = face.vertices[2] - face.vertices[1]
        diag_ind = 1 # diagonal wall index
        mirror_ind = 2 # mirror point new index
    elseif max_index == 2
        point_to_mirror = face.vertices[1]
        startpoint = face.vertices[2]
        line = face.vertices[3] - face.vertices[2]
        diag_ind = 2 # diagonal wall index
        mirror_ind = 3 # mirror point new index
    else # max_index == 3
        point_to_mirror = face.vertices[2]
        startpoint = face.vertices[3]
        line = face.vertices[1] - face.vertices[3]
        diag_ind = 3 # diagonal wall index
        mirror_ind = 4 # mirror point new index
    end
    line_midpoint = startpoint + line/2
    vec = point_to_mirror - line_midpoint
    mirrored = - vec + line_midpoint
        
    if max_index == 1
        new_points = SVector(face.vertices[1], mirrored, face.vertices[2], face.vertices[3])
        new_solid = SVector(face.solidWalls[1], face.solidWalls[1], face.solidWalls[2], face.solidWalls[3])
    elseif max_index == 2
        new_points = SVector(face.vertices[1], face.vertices[2], mirrored, face.vertices[3])
        new_solid = SVector(face.solidWalls[1], face.solidWalls[2], face.solidWalls[2], face.solidWalls[3])
    else # max_index == 3
        new_points = SVector(face.vertices[1], face.vertices[2], face.vertices[3], mirrored)
        new_solid = SVector(face.solidWalls[1], face.solidWalls[2], face.solidWalls[3], face.solidWalls[3])
    end
    # points to keep (triangle indices)
    tria_ids = setdiff([1,2,3,4], [mirror_ind]) # order of triangle points (counter clockwise) (from end of diagonal)

    # create the new face (now quadrilateral) - UPDATED with spectral support
    face2 = PolyVolume2D{G}(new_points, new_solid, n_bins, kappa_default, sigma_s_default)

    # mesh as if it was a quadrilateral
    face2 = meshQuad(face2, Ndim, Ndim)

    # calculate the projected midpoint on the diagonal
    point_vector = triangle_midPoint - startpoint # Vector from the start of the line to the point
    t = dot(point_vector, line) / dot(line, line) # Calculate the projection of point_vector onto line_vector
    t = clamp(t, 0, 1) # Clamp t to [0, 1] to ensure the nearest point is on the line segment
    nearest_point = startpoint + t * line # Calculate the nearest point

    for subface in face2.subFaces
        # check the direction of the subface, as compared to the diagonal
        costheta_subface_triMidpoint = dot(triangle_midPoint-nearest_point, subface.midPoint-nearest_point) # if > 0, same direction
        # costheta_diagonal = dot(startpoint-nearest_point, subface.midPoint-nearest_point)
        if isapprox(costheta_subface_triMidpoint, 0.0, atol=1e-6) # approx perpendicular
            # the subface is on the diagonal
            # but which side to keep?

            # find the triangle walls which are solid
            subface_solid_wall_ids = setdiff([1,2,3,4], [mirror_ind-1, mirror_ind]) # these are from the sub quadrilateral
            walls_ids = sort([subface_solid_wall_ids..., diag_ind]) # the diag is from the original triangle
            walls_solid = [i == diag_ind ? face.solidWalls[diag_ind] : subface.solidWalls[i] for i in walls_ids]

            # create the new subface - UPDATED with spectral support
            subface_points = SVector{3}(subface.vertices[tria_ids]...)
            subface_solid = SVector{3}(walls_solid...)
            subface_keeper = PolyVolume2D{G}(subface_points, subface_solid, n_bins, kappa_default, sigma_s_default)
            addSubFace!(face, subface_keeper)
        else
            # subface is not on the diagonal
            costheta_diagonal_subface = dot(subface.midPoint-nearest_point, triangle_midPoint-nearest_point)
            # costheta_diagonal_triangle = dot(startpoint-nearest_point, triangle_midPoint-nearest_point)
            if costheta_diagonal_subface > 0.0 - 1e-6
                # keep subface
                # push!(face3.subFaces, subface)
                addSubFace!(face, subface) # add the subface
            else
                # delete subface
                continue
            end
        end
    end

    return face

end

# updated 3d version of meshQuad with spectral support
function meshTriangle(face::Vector{Vector{G}}, Nx::P, Ny::P) where {G,P<:Integer}
    # save triangle midpoint
    triangle_midPoint = sum(face)/3
    
    # insert a point opposite of one the point opposite the longest edge
    norm1 = norm(face[1][1:2] - face[2][1:2])
    norm2 = norm(face[2][1:2] - face[3][1:2])
    norm3 = norm(face[3][1:2] - face[1][1:2])
    max_norm, max_index = findmax([norm1, norm2, norm3])
    
    # mirror the opposite point
    if max_index == 1
        point_to_mirror = face[3][1:2]
        startpoint = face[1][1:2]
        line = face[2][1:2] - face[1][1:2]
        diag_ind = 1 # diagonal wall index
        mirror_ind = 2 # mirror point new index
    elseif max_index == 2
        point_to_mirror = face[1][1:2]
        startpoint = face[2][1:2]
        line = face[3][1:2] - face[2][1:2]
        diag_ind = 2 # diagonal wall index
        mirror_ind = 3 # mirror point new index
    else # max_index == 3
        point_to_mirror = face[2][1:2]
        startpoint = face[3][1:2]
        line = face[1][1:2] - face[3][1:2]
        diag_ind = 3 # diagonal wall index
        mirror_ind = 4 # mirror point new index
    end
    line_midpoint = startpoint + line/2
    vec = point_to_mirror - line_midpoint
    mirrored = - vec + line_midpoint
        
    if max_index == 1
        new_points = SVector(
            Point2(face[1][1:2]...),
            Point2(mirrored...),
            Point2(face[2][1:2]...),
            Point2(face[3][1:2]...))
    elseif max_index == 2
        new_points = SVector(
            Point2(face[1][1:2]...),
            Point2(face[2][1:2]...),
            Point2(mirrored...),
            Point2(face[3][1:2]...),
        )
    else # max_index == 3
        new_points = SVector(
            Point2(face[1][1:2]...),
            Point2(face[2][1:2]...),
            Point2(face[3][1:2]...),
            Point2(mirrored...),
        )
    end
    # points to keep (triangle indices)
    tria_ids = setdiff([1,2,3,4], [mirror_ind]) # order of triangle points (counter clockwise) (from end of diagonal)

    # create the new face (now quadrilateral)
    new_solid = SVector(true, true, true, true)
    face2 = PolyVolume2D{G}(new_points, new_solid, 1, 1.0, 1.0)

    # mesh as if it was a quadrilateral
    face2 = meshQuad(face2, Nx, Ny)

    # calculate the projected midpoint on the diagonal
    point_vector = triangle_midPoint[1:2] - startpoint # Vector from the start of the line to the point
    t = dot(point_vector, line) / dot(line, line) # Calculate the projection of point_vector onto line_vector
    t = clamp(t, 0, 1) # Clamp t to [0, 1] to ensure the nearest point is on the line segment
    nearest_point = startpoint + t * line # Calculate the nearest point

    # Create the new points
    point1 = Point3[]
    point2 = Point3[]
    point3 = Point3[]
    point4 = Point3[]
    for subface in face2.subFaces
        # check the direction of the subface, as compared to the diagonal
        costheta_subface_triMidpoint = dot(triangle_midPoint[1:2]-nearest_point, subface.midPoint-nearest_point) # if > 0, same direction
        # costheta_diagonal = dot(startpoint-nearest_point, subface.midPoint-nearest_point)
        if isapprox(costheta_subface_triMidpoint, 0.0, atol=1e-6) # approx perpendicular
            # the subface is on the diagonal
            # but which side to keep?

            # find the triangle walls which are solid
            subface_solid_wall_ids = setdiff([1,2,3,4], [mirror_ind-1, mirror_ind]) # these are from the sub quadrilateral
            walls_ids = sort([subface_solid_wall_ids..., diag_ind]) # the diag is from the original triangle
            walls_solid = [true, true, true] # [i == diag_ind ? face.solidWalls[diag_ind] : subface.solidWalls[i] for i in walls_ids]

            # create the new subface
            subface_points = SVector{3}(subface.vertices[tria_ids]...)
            subface_solid = SVector{3}(walls_solid...)
            subface_keeper = PolyVolume2D{G}(subface_points, subface_solid, 1, 1.0, 1.0)
            push!(point1, Point3(subface_keeper.vertices[1]..., 0.0))
            push!(point2, Point3(subface_keeper.vertices[2]..., 0.0))
            push!(point3, Point3(subface_keeper.vertices[3]..., 0.0))
            push!(point4, Point3(subface_keeper.vertices[3]..., 0.0))
        else
            # subface is not on the diagonal
            costheta_diagonal_subface = dot(subface.midPoint-nearest_point, triangle_midPoint[1:2]-nearest_point)
            # costheta_diagonal_triangle = dot(startpoint-nearest_point, triangle_midPoint-nearest_point)
            if costheta_diagonal_subface > 0.0 - 1e-6
                # keep subface
                push!(point1, Point3(subface.vertices[1]..., 0.0))
                push!(point2, Point3(subface.vertices[2]..., 0.0))
                push!(point3, Point3(subface.vertices[3]..., 0.0))
                push!(point4, Point3(subface.vertices[4]..., 0.0))
            else
                # delete subface
                continue
            end
        end
    end

    return point1, point2, point3, point4

end