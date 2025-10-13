# Updated meshQuad with spectral support
function meshQuad(face::PolyFace2D{G}, Nx::P, Ny::P) where {G,P<:Integer}
    # Determine spectral mode from parent face
    is_spectral = isa(face.kappa_g, Vector)
    n_bins = is_spectral ? length(face.kappa_g) : 1
    
    # Get default extinction values from parent face
    kappa_default = is_spectral ? face.kappa_g[1] : face.kappa_g
    sigma_s_default = is_spectral ? face.sigma_s_g[1] : face.sigma_s_g
    
    # Define points in enclosure.
    Apoint = face.vertices[1]
    Bpoint = face.vertices[2]
    Cpoint = face.vertices[3]
    Dpoint = face.vertices[4]
    # x- and y-coordinates of this enclosure.
    xs = SVector(Apoint[1], Bpoint[1], Cpoint[1], Dpoint[1], Apoint[1])
    ys = SVector(Apoint[2], Bpoint[2], Cpoint[2], Dpoint[2], Apoint[2])

    # Split each sub-enclosure into a specified number cells.

    # Initialize arrays to hold information about each cell.
    xPoints = zeros(G, Nx+1,Ny+1)
    yPoints = zeros(G, Nx+1,Ny+1)
    point1 = Matrix{Point2{G}}(undef, Nx, Ny)
    point2 = Matrix{Point2{G}}(undef, Nx, Ny)
    point3 = Matrix{Point2{G}}(undef, Nx, Ny)
    point4 = Matrix{Point2{G}}(undef, Nx, Ny)

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

        # reset the solid walls vector
        solid_walls = [false, false, false, false]

        # check for solid walls to remain solid
        if m == 1
            solid_walls[1] = face.solidWalls[1]
        elseif m == Ny
            solid_walls[3] = face.solidWalls[3]
        end

        for n = 1:Nx # loop over x

            # reset the solid walls vector
            solid_walls = [solid_walls[1], false, solid_walls[3], false]

            # check for solid walls to remain solid
            if n == 1
                solid_walls[4] = face.solidWalls[4]
            elseif n == Nx
                solid_walls[2] = face.solidWalls[2]
            end

            # convert to static array
            solid_walls_in = SVector(solid_walls[1], solid_walls[2], solid_walls[3], solid_walls[4])

            # Create a new PolyFace2D subface - UPDATED with spectral support
            point1 = Point2(xPoints[n,m], yPoints[n,m])
            point2 = Point2(xPoints[n+1,m], yPoints[n+1,m])
            point3 = Point2(xPoints[n+1,m+1], yPoints[n+1,m+1])
            point4 = Point2(xPoints[n,m+1], yPoints[n,m+1])
            subface = PolyFace2D{G}(SVector(point1, point2, point3, point4), solid_walls_in, 
                                     n_bins, kappa_default, sigma_s_default)

            # insert it
            # push!(face.subFaces, subface)
            addSubFace!(face, subface)

        end
    end

    return face
end

# Updated meshTriangle with spectral support
function meshTriangle(face::PolyFace2D{G}, Ndim::P) where {G, P<:Integer}
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
    face2 = PolyFace2D{G}(new_points, new_solid, n_bins, kappa_default, sigma_s_default)

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
            subface_keeper = PolyFace2D{G}(subface_points, subface_solid, n_bins, kappa_default, sigma_s_default)
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

function addSubFace!(superFace::PolyFace2D{G}, subFace::PolyFace2D{G}) where {G}

    # next, update the properties of the added subface
    # the subface inherits the properties of the superface

    # inherit extinction properties - direct assignment/copy
    if isa(superFace.kappa_g, Vector)
        subFace.kappa_g .= superFace.kappa_g
        subFace.sigma_s_g .= superFace.sigma_s_g
    else
        subFace.kappa_g = superFace.kappa_g
        subFace.sigma_s_g = superFace.sigma_s_g
    end

    # inherit volume state variables - direct assignment/copy
    inherit_volume_property!(superFace, subFace, :j_g)
    inherit_volume_property!(superFace, subFace, :g_a_g)
    inherit_volume_property!(superFace, subFace, :e_g)
    inherit_volume_property!(superFace, subFace, :r_g)
    inherit_volume_property!(superFace, subFace, :g_g)
    inherit_volume_property!(superFace, subFace, :i_g)
    inherit_volume_property!(superFace, subFace, :q_in_g)
    inherit_volume_property!(superFace, subFace, :q_g)
    inherit_volume_property!(superFace, subFace, :T_in_g)
    inherit_volume_property!(superFace, subFace, :T_g)

    # then for walls
    for i in 1:length(subFace.vertices) # loop over walls of subface (same as number of vertices)
        if (length(superFace.vertices) == 3 && length(subFace.vertices) == 4) && subFace.solidWalls[i] == true
            # if superFace is a triangle and subFace is a quadrilateral
            # inherits the properties of the superface
            
            # two walls of the quadrilateral need to inherit from the diagonal of the triangle

            # find the diagonal (should still work for equilateral triangles)
            norm1 = norm(superFace.vertices[1]-superFace.vertices[2])
            norm2 = norm(superFace.vertices[2]-superFace.vertices[3])
            norm3 = norm(superFace.vertices[3]-superFace.vertices[1])
            diag_len, diag_index = findmax([norm1, norm2, norm3])
            diag_normal = superFace.outwardNormals[diag_index]
            
            # check if the outwardNormal of the subface is in the same direction as the diagonal outwardNormal
            diag_dirs = [dot(subFace.outwardNormals[j], diag_normal) for j in 1:4] # positive for sides which inherit from the diagonal
            # and that the diagonal is solid
            # need to +1 everywhere to account for i=1 is volume
            if diag_dirs[i] > 0.0 && superFace.solidWalls[diag_index] == true
                # inherits the properties of the superface (the diagonal)
                inheritSurfaceProperties!(superFace, subFace; from=diag_index, to=i)

            elseif diag_dirs[i] < 0.0 && superFace.solidWalls[diag_index] == true
                # if the quadrilateral subface wall is solid, but should not inherit from the superFace triangle diagonal
                # which index to inherit from? There are four cases
                # 1) subface sides 1 and 2 inherits from diagonal, subface side 3 and 4 inherits from triangle sides 2 and 3, valid
                if diag_dirs[1] > 0.0 && diag_dirs[2] > 0.0
                    # subface side 3 inherits from triangle side 2
                    inheritSurfaceProperties!(superFace, subFace; from=2, to=3)

                    # subface side 4 inherits from triangle side 3
                    inheritSurfaceProperties!(superFace, subFace; from=3, to=4)
                end

                # 2) subface sides 2 and 3 inherits from diagonal, subface side 1 and 4 inherits from triangle sides 1 and 3, valid
                if diag_dirs[2] > 0.0 && diag_dirs[3] > 0.0
                    # subface side 1 inherits from triangle side 1
                    inheritSurfaceProperties!(superFace, subFace; from=1, to=1)

                    # subface side 4 inherits from triangle side 3
                    inheritSurfaceProperties!(superFace, subFace; from=3, to=4)
                end

                # 3) subface sides 3 and 4 inherits from diagonal, subface side 1 and 2 inherits from triangle sides 1 and 2, valid
                if diag_dirs[3] > 0.0 && diag_dirs[4] > 0.0
                    # subface side 1 inherits from triangle side 1
                    inheritSurfaceProperties!(superFace, subFace; from=1, to=1)

                    # subface side 2 inherits from triangle side 2
                    inheritSurfaceProperties!(superFace, subFace; from=2, to=2)                    
                end

                # 4) subface sides 1 and 4 inherits from diagonal, subface side 2 and 3 inherits from triangle sides 1 and 2, valid
                if diag_dirs[4] > 0.0 && diag_dirs[1] > 0.0
                    # subface side 2 inherits from triangle side 1
                    inheritSurfaceProperties!(superFace, subFace; from=1, to=2)

                    # subface side 4 inherits from triangle side 3
                    inheritSurfaceProperties!(superFace, subFace; from=2, to=3)
                end

            end

        elseif (length(superFace.vertices) == length(subFace.vertices)) && subFace.solidWalls[i] == true
            # if both are quadrilaterals or both are triangles
            inheritSurfaceProperties!(superFace, subFace; from=i, to=i)            
        else
            # wall is not solid, no need to inherit
            continue
        end
    end

    # push the subFace to the list of subFaces
    push!(superFace.subFaces, subFace)

end

# Helper function to inherit volume properties - direct copy
function inherit_volume_property!(superFace::PolyFace2D{G}, subFace::PolyFace2D{G}, property::Symbol) where {G}
    if property == :q_in_g || property == :q_g
        super_val = getfield(superFace, property)
        sub_val = getfield(subFace, property)

        # source flux should reduce by area ratio
        sub_val = super_val*subFace.volume/superFace.volume
        setfield!(subFace, property, sub_val)
    else
        super_val = getfield(superFace, property)
        sub_val = getfield(subFace, property)
    
        if isa(super_val, Vector)
            # Spectral - copy entire vector
            sub_val .= super_val
        else
            # Grey - direct assignment
            setfield!(subFace, property, super_val)
        end
    end
end

# Updated inheritSurfaceProperties! function with direct inheritance
function inheritSurfaceProperties!(superFace::PolyFace2D{G}, subFace::PolyFace2D{G}; from=i, to=j) where {G}
    # the subFace inherits the properties of the superFace
    
    # boundary properties (unchanged - always scalar)
    subFace.epsilon[to] = superFace.epsilon[from]
    
    # wall state variables - direct inheritance
    inherit_wall_property!(superFace, subFace, :j_w, from, to)
    inherit_wall_property!(superFace, subFace, :g_a_w, from, to)
    inherit_wall_property!(superFace, subFace, :e_w, from, to)
    inherit_wall_property!(superFace, subFace, :r_w, from, to)
    inherit_wall_property!(superFace, subFace, :g_w, from, to)
    inherit_wall_property!(superFace, subFace, :i_w, from, to)
    inherit_wall_property!(superFace, subFace, :q_in_w, from, to)
    inherit_wall_property!(superFace, subFace, :q_w, from, to)
    inherit_wall_property!(superFace, subFace, :T_in_w, from, to)
    inherit_wall_property!(superFace, subFace, :T_w, from, to)
end

# Helper function to inherit wall properties - direct copy
function inherit_wall_property!(superFace::PolyFace2D{G}, subFace::PolyFace2D{G}, 
                               property::Symbol, from::Int, to::Int) where {G}
    if property == :q_in_w || property == :q_w
        
        super_wall_array = getfield(superFace, property)
        sub_wall_array = getfield(subFace, property)
        
        super_val = super_wall_array[from]
        sub_val = sub_wall_array[to]
        
        # source flux should reduce by area ratio
        sub_wall_array[to] = super_val*subFace.area[to]/superFace.area[from]
    else
        super_wall_array = getfield(superFace, property)
        sub_wall_array = getfield(subFace, property)
        
        super_val = super_wall_array[from]
        sub_val = sub_wall_array[to]
        
        if isa(super_val, Vector)
            # Spectral - copy entire vector
            sub_val .= super_val
        else
            # Grey - direct assignment
            sub_wall_array[to] = super_val
        end
    end
end