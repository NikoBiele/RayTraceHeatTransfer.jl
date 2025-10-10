function isotropicScatter(nudge::G) where {G}
    theta = acos(2*rand() - 1)
    phi = 2π * rand()
    return Point2{G}(G(sin(theta)*cos(phi)), G(cos(theta)))
end