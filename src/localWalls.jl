function localWalls(mesh::TracingMesh, N_subs_count::Int64)
    # this function determines a point on each wall of the current cell
    # as well as normal vectors of each wall
    # this information fully describe the walls of the current cell
    
    # calculate points on each wall    
    pointOne = mesh.point1_coarse[1,N_subs_count]
    pointTwo = mesh.point2_coarse[1,N_subs_count]
    pointThree = mesh.point3_coarse[1,N_subs_count]
    pointFour = mesh.point4_coarse[1,N_subs_count]
    
    # for output readability
    wallPointBottom = pointOne
    wallPointRight = pointTwo
    wallPointTop = pointThree
    wallPointLeft = pointFour

    # now calculate normal vectors for a complete description
    bottomWallNormal = SVector(0.0, 1.0) # top and bottom always vertical
    topWallNormal = SVector(0.0, -1.0)
    
    # the right-hand wall
    firstPoint = pointTwo
    secondPoint = pointThree
    if abs(firstPoint[1]-secondPoint[1]) < 1e-3 # if the change in x is less than one mm then we say it's a vertical line

        # normal vector
        rightWallNormal = SVector(-1.0, 0.0)

    else # the line is not vertical so we construct the line

        a = (secondPoint[2]-firstPoint[2])/(secondPoint[1]-firstPoint[1])
        dx = abs(secondPoint[1]-firstPoint[1])
        dy = abs(secondPoint[2]-firstPoint[2])
        if a < 0 # if line is decreasing
            rightWallNormal = SVector{2}([-dy, -dx]./norm(secondPoint-firstPoint)) # works
        elseif a > 0  # if line is increasing
            rightWallNormal = SVector{2}([-dy, dx]./norm(secondPoint-firstPoint)) # works
        end

    end

    # left-hand wall
    firstPoint = pointFour
    secondPoint = pointOne
    if abs(firstPoint[1]-secondPoint[1]) < 1e-3 # if the change in x is less than one mm then we say it's a vertical line

        leftWallNormal = SVector(1.0, 0.0)              

    else # the line is not vertical so we construct the line

        a = (secondPoint[2]-firstPoint[2])/(secondPoint[1]-firstPoint[1])
        dx = abs(secondPoint[1]-firstPoint[1])
        dy = abs(secondPoint[2]-firstPoint[2])
        if a < 0 # if line is decreasing
            leftWallNormal = SVector{2}([dy, dx]./norm(secondPoint-firstPoint)) # works
        elseif a > 0  # if line is increasing
            leftWallNormal = SVector{2}([dy, -dx]./norm(secondPoint-firstPoint)) # works
        end

    end

    return wallPointBottom, wallPointRight, wallPointTop, wallPointLeft, bottomWallNormal, rightWallNormal, topWallNormal, leftWallNormal
end