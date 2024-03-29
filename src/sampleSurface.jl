function sampleSurface(mesh::TracingMesh, xCountSample::Int64, yCountSample::Int64,
                        sampleLeftRight::Bool, sampleTopBottom::Bool)

    # This function samples an emission position on a 2D surface (a line) as well as an emission direction,
    # given the lines endpoints and position in the mesh.
    # The function samples from a diffuse (Lambertian) distribution and rotates the sample
    # from local to global coordinates using a rotation matrix.
    # This function was written with a separate
    # method for each combination to make it easier to debug and validate.

    # Set (m,n), the coordinates for the cell being sampled.
    m = xCountSample
    n = yCountSample

    # From information on position in the mesh, identify solid walls bounding this cell.
    if xCountSample == 1
        solidWallx = 4
    elseif xCountSample == mesh.Nx
        solidWallx = 2
    else
        solidWallx = 0 # error protection
    end
    if yCountSample == 1
        solidWally = 1
    elseif yCountSample == mesh.N_subs*mesh.Ny
        solidWally = 3
    else
        solidWally = 0 # error protection
    end

    # Unit vectors in global coordinate system.
    xVecGlobal = SVector(1.0, 0.0)
    yVecGlobal = SVector(0.0, 1.0) # Used for rotating from local to global with rotation matrix.

    # If we sample the top and bottom we assume that the wall is horisontal.
    # Then we take a Lambertian (diffuse) sample and we do not use a rotation matrix,
    # we only flip the direction for the top wall (by taking the negative vector).
    if sampleTopBottom == true
        if solidWally == 1 # bottom wall

            # get bounding points
            firstPoint = SVector{2}(mesh.point1[m,n])
            secondPoint = SVector{2}(mesh.point2[m,n])

            # Construct line and sample an emission point and direction.
            x0 = minimum([firstPoint[1], secondPoint[1]])
            R_posx = rand()
            dx = abs(firstPoint[1]-secondPoint[1])
            x = x0 + dx*R_posx
            y = firstPoint[2]
            point = SVector(x, y+1e-9) # sampled position (slightly above wall)
            i1_loc = lambertSample3D() # sample direction of ray
            i1 = SVector{2}(i1_loc)

        elseif solidWally == 3 # top wall

            # get bounding points
            firstPoint = SVector{2}(mesh.point3[m,n])
            secondPoint = SVector{2}(mesh.point4[m,n])

            # construct line and sample an emission point and direction
            R_posx = rand()
            x0 = minimum([firstPoint[1], secondPoint[1]])
            dx = abs(firstPoint[1]-secondPoint[1])
            x = x0 + dx*R_posx
            y = firstPoint[2]
            point = SVector(x, y-1e-9) # sampled position (slightly below wall)
            i1_loc = lambertSample3D() # sample direction of ray
            i1 = SVector{2}(-i1_loc)

        end
    end
    
    # If we sample left or right wall we need to check for
    # vertical first (to avoid division by zero)
    if sampleLeftRight == true
        if solidWallx == 2 # right wall
            
            # get bounding points
            firstPoint = SVector{2}(mesh.point2[m,n])
            secondPoint = SVector{2}(mesh.point3[m,n])

            # If the change in x is less than one mm then we say it's a vertical line
            if abs(firstPoint[1]-secondPoint[1]) < 1e-3

                # then we only have to sample the y-position
                dy = abs(firstPoint[2]-secondPoint[2])
                R_posy = rand() # sample for position
                y = firstPoint[2] + dy*R_posy
                x = minimum([firstPoint[1], secondPoint[1]]) # must be minimum since we're at wall right (to ensure we're inside domain)
                point = SVector(x-1e-9, y)
                # sample direction of ray (local coordinate system)
                i1_loc = lambertSample3D()
                # unit vectors in local coordinate system
                xVecLocal = SVector{2}([0.0, 1.0])
                yVecLocal = SVector{2}([-1.0, 0.0])
                # rotate according to rotation matrix 
                RotationMatrix = SMatrix{2,2}([dot(xVecGlobal, xVecLocal) dot(xVecGlobal, yVecLocal);
                                                dot(yVecGlobal, xVecLocal) dot(yVecGlobal, yVecLocal)])
                i1 = SVector{2}(RotationMatrix*i1_loc) # emission direction (global coordinate system)

            else # the line is not vertical so we construct the line and sample it

                # right bounding line
                # first we find the point from which we emit
                a = (secondPoint[2] - firstPoint[2])/(secondPoint[1] - firstPoint[1])
                b = firstPoint[2] - a*firstPoint[1]
                # sample the x position and get the y position
                R_posx = rand()
                x0 = minimum([firstPoint[1], secondPoint[1]]) # must be minimum to construct line
                dx = abs(secondPoint[1]-firstPoint[1])
                x = x0 + dx*R_posx
                y = a*x + b
                point = SVector(x-1e-9, y)

                # Sample direction of ray, first get the local vectors of the wall
                # using this validated method (one method for each combination for debugging)
                dy = abs(secondPoint[2]-firstPoint[2])
                if a < 0 # if line is decreasing
                    xVecLocal = SVector{2}([-dx, dy]./norm(secondPoint-firstPoint)) # works
                    yVecLocal = SVector{2}([-dy, -dx]./norm(secondPoint-firstPoint)) # works
                elseif a > 0  # if line is increasing
                    xVecLocal = SVector{2}([dx, dy]./norm(secondPoint-firstPoint)) # works
                    yVecLocal = SVector{2}([-dy, dx]./norm(secondPoint-firstPoint)) # works
                end
                # sample direction of ray (local coordinate system)
                i1_loc = lambertSample3D()
                # rotate according to rotation matrix 
                RotationMatrix = SMatrix{2,2}([dot(xVecGlobal, xVecLocal) dot(xVecGlobal, yVecLocal);
                                                dot(yVecGlobal, xVecLocal) dot(yVecGlobal, yVecLocal)])
                i1 = SVector{2}(RotationMatrix*i1_loc)
            end

        elseif solidWallx == 4 # then sample the left wall

            # get bounding points
            firstPoint = SVector{2}(mesh.point4[m,n])
            secondPoint = SVector{2}(mesh.point1[m,n])

            if abs(firstPoint[1]-secondPoint[1]) < 1e-3 # if the change in x is less than one mm then we say it's a vertical line

                # then we only have to sample the y-position
                dy = abs(firstPoint[2]-secondPoint[2])
                R_posy = rand() # sample for position
                y = secondPoint[2] + dy*R_posy
                x = maximum([firstPoint[1], secondPoint[1]]) # must be maximum since we're at wall left (to ensure we're inside domain)
                point = SVector(x+1e-9, y)
                # sample direction of ray (local coordinate system)
                i1_loc = lambertSample3D()
                # unit vectors in local coordinate system
                xVecLocal = SVector{2}([0.0, -1.0])
                yVecLocal = SVector{2}([1.0, 0.0])
                # rotate according to rotation matrix 
                RotationMatrix = SMatrix{2,2}([dot(xVecGlobal, xVecLocal) dot(xVecGlobal, yVecLocal);
                                                dot(yVecGlobal, xVecLocal) dot(yVecGlobal, yVecLocal)])
                i1 = SVector{2}(RotationMatrix*i1_loc)

            else # the line is not vertical so we construct the line and sample it

                # first we find the point from which we emit
                a = (secondPoint[2]-firstPoint[2])/(secondPoint[1]-firstPoint[1])
                b = firstPoint[2] - a * firstPoint[1]
                # sample the x position and get the y position
                R_posx = rand()
                x0 = minimum([firstPoint[1], secondPoint[1]]) # must be minimum to construct line
                dx = abs(secondPoint[1]-firstPoint[1])
                x = x0 + dx*R_posx
                y = a*x + b
                point = SVector(x+1e-9, y)

                # Sample direction of ray, first get the local vectors of the wall
                # using this validated method (one method for each combination for debugging)
                dy = abs(secondPoint[2]-firstPoint[2])
                if a < 0 # if line is decreasing
                    xVecLocal = [dx, -dy]./norm(secondPoint-firstPoint) # works
                    yVecLocal = [dy, dx]./norm(secondPoint-firstPoint) # works
                elseif a > 0  # if line is increasing
                    xVecLocal = [-dx, -dy]./norm(secondPoint-firstPoint) # works
                    yVecLocal = [dy, -dx]./norm(secondPoint-firstPoint) # works
                end
                # sample direction of ray (local coordinate system)
                i1_loc = lambertSample3D()
                # rotate according to rotation matrix
                RotationMatrix = SMatrix{2,2}([dot(xVecGlobal, xVecLocal) dot(xVecGlobal, yVecLocal);
                                                dot(yVecGlobal, xVecLocal) dot(yVecGlobal, yVecLocal)])
                i1 = SVector{2}(RotationMatrix*i1_loc)

            end
        
        end
    
    end

    return point, i1
    
end