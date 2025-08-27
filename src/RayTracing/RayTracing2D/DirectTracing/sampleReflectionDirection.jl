function sample_reflection_direction(normal::Point2{Float64})
    # Sample direction in local coordinate system
    i1_loc = lambertSample3D()
    
    # Create local coordinate system
    xVecLocal = SVector{2}([normal[2], -normal[1]])
    yVecLocal = normal
    
    # Rotation matrix from local to global coordinate system
    RotationMatrix = SMatrix{2,2}([dot(xVecGlobal, xVecLocal) dot(xVecGlobal, yVecLocal);
                                   dot(yVecGlobal, xVecLocal) dot(yVecGlobal, yVecLocal)])
    
    # Transform to global coordinate system
    i1_init = SVector{2}(RotationMatrix * i1_loc)
    
    return Point2(i1_init[1], i1_init[2])
end