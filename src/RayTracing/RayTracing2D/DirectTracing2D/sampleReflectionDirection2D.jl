function sampleReflectionDirection2D(normal::Point2{G}, rng::AbstractRNG) where {G}
    # sample the direction in the wall's local coordinate system
    # (x along the wall, y along the inward normal)
    i1_loc = lambertSample2D(rng, G)

    # local unit vectors in global coordinates: x = (n2, -n1), y = n. The global basis is the
    # identity, so the rotation matrix has these as its columns; the product is written out
    return Point2{G}( normal[2]*i1_loc[1] + normal[1]*i1_loc[2],
                     -normal[1]*i1_loc[1] + normal[2]*i1_loc[2])
end