function sample_reflection_direction(normal::Point2{Float64})
    # Sample direction in local coordinate system
    i1_loc = lambertSample3D()
    
    # Create local coordinate system
    xVecLocal2D = SVector{2}([normal[2], -normal[1]])
    yVecLocal2D = normal
    
    # Rotation matrix from local to global coordinate system
    RotationMatrix2D = SMatrix{2,2}([dot(xVecGlobal2D, xVecLocal2D) dot(xVecGlobal2D, yVecLocal2D);
                                   dot(yVecGlobal2D, xVecLocal2D) dot(yVecGlobal2D, yVecLocal2D)])
    
    # Transform to global coordinate system
    i1_init = SVector{2}(RotationMatrix2D * i1_loc)
    
    return Point2(i1_init[1], i1_init[2])
end