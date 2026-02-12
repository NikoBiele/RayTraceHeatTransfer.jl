# updated 3d view factor version of meshQuad with spectral support
function meshQuad(face::Vector{Vector{G}}, Nx::P, Ny::P) where {G,P<:Integer}
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

# Updated 2d view factor version of meshQuad with spectral support
function meshQuad(volume::PolyVolume2D{G}, Nx::P, Ny::P) where {G,P<:Integer}
    # Determine spectral mode from parent volume
    is_spectral = isa(volume.kappa_g, Vector)
    n_bins = is_spectral ? length(volume.kappa_g) : 1
    
    # Get default extinction values from parent volume
    kappa_default = is_spectral ? volume.kappa_g[1] : volume.kappa_g
    sigma_s_default = is_spectral ? volume.sigma_s_g[1] : volume.sigma_s_g
    
    # Define points in enclosure.
    Apoint = volume.vertices[1]
    Bpoint = volume.vertices[2]
    Cpoint = volume.vertices[3]
    Dpoint = volume.vertices[4]
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
            solid_walls[1] = volume.solidWalls[1]
        elseif m == Ny
            solid_walls[3] = volume.solidWalls[3]
        end

        for n = 1:Nx # loop over x

            # reset the solid walls vector
            solid_walls = [solid_walls[1], false, solid_walls[3], false]

            # check for solid walls to remain solid
            if n == 1
                solid_walls[4] = volume.solidWalls[4]
            elseif n == Nx
                solid_walls[2] = volume.solidWalls[2]
            end

            # convert to static array
            solid_walls_in = SVector(solid_walls[1], solid_walls[2], solid_walls[3], solid_walls[4])

            # Create a new PolyVolume2D subvolume - UPDATED with spectral support
            point1 = Point2(xPoints[n,m], yPoints[n,m])
            point2 = Point2(xPoints[n+1,m], yPoints[n+1,m])
            point3 = Point2(xPoints[n+1,m+1], yPoints[n+1,m+1])
            point4 = Point2(xPoints[n,m+1], yPoints[n,m+1])
            subvolume = PolyVolume2D{G}(SVector(point1, point2, point3, point4), solid_walls_in, 
                                     n_bins, kappa_default, sigma_s_default)

            # insert it
            # push!(volume.subFaces, subvolume)
            addSubVolume!(volume, subvolume)

        end
    end

    return volume
end
