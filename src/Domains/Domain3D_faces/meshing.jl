# Mesh a quadrilateral
function meshQuad(face::Vector{Vector{G}}, Nx::Int, Ny::Int) where {G<:AbstractFloat}
    faces = Vector{Vector{G}}[] # for output

    # Define points in enclosure.
    Apoint = face[1]
    Bpoint = face[2]
    Cpoint = face[3]
    Dpoint = face[4]
    # x- and y-coordinates of this enclosure.
    xs = SVector(Apoint[1], Bpoint[1], Cpoint[1], Dpoint[1], Apoint[1])
    ys = SVector(Apoint[2], Bpoint[2], Cpoint[2], Dpoint[2], Apoint[2])

    # Split each sub-enclosure into a specified number cells.

    # Initialize arrays to hold information about each cell.
    xPoints = zeros(G, Nx+1, Ny+1)
    yPoints = zeros(G, Nx+1, Ny+1)
    point1 = Point3[] # Matrix{Point3{G}}(undef, Nx, Ny)
    point2 = Point3[] # Matrix{Point3{G}}(undef, Nx, Ny)
    point3 = Point3[] # Matrix{Point3{G}}(undef, Nx, Ny)
    point4 = Point3[] # Matrix{Point3{G}}(undef, Nx, Ny)

    # Improved algorithm!

    # Get change of coordinates at outer dimensions.
    # x first.
    deltaXbot = xs[2]-xs[1]
    deltaXtop = xs[4]-xs[3]
    deltaXleft = xs[5]-xs[4]
    # then y.
    deltaYBot = ys[1]-ys[2]
    deltaYRight = ys[2] - ys[3]
    deltaYLeft = ys[4] - ys[1]

    # Loop to create sub-divisions (cells) of each sub-enclosure.
    for m = 1:Ny+1 # loop over y

        # We move the reference point in each iteration.
        refmoveXleft = (m-1)*deltaXleft/Ny
        refmoveXright = deltaXbot - (m-1)*(deltaXbot+deltaXtop)/Ny

        # Start at the bottom left point in aNdim enclosure.
        for n = 1:Nx+1 # loop over x

            # Change in x and y.
            # Add or subtract a bit of y.
            refmoveYdown = (n-1)*deltaYBot/Nx
            refmoveYup = deltaYLeft - (n-1)*(deltaYLeft+deltaYRight)/Nx

            # First term = bottom left point.
            # Second term = how much we move the reference.
            # Third term = how much we move right in each iteration.
            xPoints[n,m] = xs[1] - refmoveXleft + (n-1)*(refmoveXright)/Nx
            yPoints[n,m] = ys[1] - refmoveYdown + (m-1)*(refmoveYup)/Ny
        end
    end
    
    # Put the points into the matrices to work with ray tracing.
    for m = 1:Ny # loop over y
        for n = 1:Nx # loop over x

            # Create the new points
            push!(point1, Point3(xPoints[n,m], yPoints[n,m], 0.0))
            push!(point2, Point3(xPoints[n+1,m], yPoints[n+1,m], 0.0))
            push!(point3, Point3(xPoints[n+1,m+1], yPoints[n+1,m+1], 0.0))
            push!(point4, Point3(xPoints[n,m+1], yPoints[n,m+1], 0.0))
        end
    end

    return point1, point2, point3, point4
end

function meshTriangle(face::Vector{Vector{G}}, Nx::Int, Ny::Int) where {G<:AbstractFloat}
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
    face2 = PolyFace2D{G,G}(new_points, new_solid)

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
            subface_keeper = PolyFace2D{G,G}(subface_points, subface_solid)
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

function mesh3D_faces(coords, faces, Ndim)

    num_faces = size(faces, 1)
    points_per_face = size(faces, 2)

    # Project vertices to xy-Plane
    projected_faces = [project_face_flat(coords, faces[i, :]) for i in 1:num_faces]

    # mesh each face in the xy-plane
    mesh_pre_project = [points_per_face == 4 ? meshQuad(projected_faces[i][1],Ndim,Ndim) :
                        meshTriangle(projected_faces[i][1],Ndim,Ndim) for i in 1:num_faces]

    # then project each face back
    mesh_post_project = [project_face_back(mesh_pre_project[i], projected_faces[i][2], projected_faces[i][3]) for i in 1:num_faces]

    return mesh_post_project
end