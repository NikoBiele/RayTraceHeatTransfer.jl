function isotropicScatter()
    theta = acos(2rand() - 1)
    phi = 2π * rand()
    return Point2(sin(theta)*cos(phi), cos(theta))
end