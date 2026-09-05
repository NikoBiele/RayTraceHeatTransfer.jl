function isotropicScatter2D(nudge::G, rng::AbstractRNG) where {G}
    theta = acos(2*rand(rng, G) - 1)
    phi = 2π * rand(rng, G)
    return Point2{G}(sin(theta)*cos(phi), cos(theta))
end