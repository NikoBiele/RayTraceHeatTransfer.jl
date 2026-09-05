function lambertSample2D(rng::AbstractRNG, G)
    R_angle1 = rand(rng, G)
    cosTheta = sqrt(R_angle1)
    sinTheta = sqrt(1.0 - cosTheta^2)
    psi = 2*G(pi)*rand(rng, G)
    
    xdir = sinTheta*cos(psi)
    zdir = cosTheta
    
    return Point2(xdir, zdir)
end