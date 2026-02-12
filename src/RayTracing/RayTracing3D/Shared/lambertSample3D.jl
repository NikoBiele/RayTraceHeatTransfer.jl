function lambertSample3D()
    R_angle1 = Float32(rand())
    cosTheta = sqrt(R_angle1)
    sinTheta = sqrt(1.0 - cosTheta^2)
    psi = 2*pi*Float32(rand())
    
    xdir = sinTheta*cos(psi)
    zdir = cosTheta
    
    return Point2(xdir, zdir)
end