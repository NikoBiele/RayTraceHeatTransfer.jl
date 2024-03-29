function lambertSample3D()
    # This function generates a local diffuse (Lambertian) emission sample from a wall.
    # The generated sample is 3D but is converted to 2D
    # by projecting it onto the 2D plane (by setting third component to zero).
    # This sample correspond to the local coordinates of the given wall.
    # It is later rotated to global coordinates.

    R_angle1 = rand() # sample for finding angle
    sinTheta = sqrt(R_angle1)
    cosTheta = sqrt(1.0-(sinTheta^2))
    R_angle2 = rand()
    psi = 2*pi*R_angle2
    zdir = cosTheta # vertical direction (local)
    xdir = sinTheta*cos(psi) # horisontal direction (local)
    i1 = SVector(xdir, zdir)

    return i1
end
