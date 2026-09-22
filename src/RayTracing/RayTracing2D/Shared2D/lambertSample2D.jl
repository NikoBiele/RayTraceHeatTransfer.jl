function lambertSample2D(rng::AbstractRNG, G)
    R_angle1 = rand(rng, G)
    cosTheta = sqrt(R_angle1)
    sinTheta = sqrt(one(G) - R_angle1)      # = sqrt(1 - cosTheta^2), kept in the mesh's float type
    psi = 2*G(pi)*rand(rng, G)
    
    xdir = sinTheta*cos(psi)
    zdir = cosTheta
    
    return Point2{G}(xdir, zdir)
end