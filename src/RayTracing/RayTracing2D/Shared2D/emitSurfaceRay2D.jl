function emitSurfaceRay2D(face::PolyVolume2D{G}, wall_index::P, nudge::G,
                        rng::AbstractRNG) where {G, P<:Integer}
    # println("Emitting surface ray, type of G is $(typeof(G))")

    p1, p2 = face.vertices[wall_index], face.vertices[mod1(wall_index+1, length(face.vertices))]
    R = rand(rng, G) # Convert random number to the appropriate type G

    # sample position to find emission point
    p = p1 + (p2 - p1) * R
    # nudge the point a tiny bit towards the midpoint, to ensure we are inside cell
    p = p + (face.midPoint - p) * nudge
    p = Point2{G}(p[1], p[2])

    # sample direction of ray (local coordinate system)
    i1_loc = lambertSample2D(rng, G) # Convert the result to type G if needed
    
    # unit vectors in local coordinate system
    xVecLocal = normalize(p2-p1)
    yVecLocal = Point2{G}(-xVecLocal[2], xVecLocal[1])

    # rotate to the global coordinate system. The global basis is the identity:
    i1 = Point2{G}(xVecLocal[1]*i1_loc[1] + yVecLocal[1]*i1_loc[2],
                   xVecLocal[2]*i1_loc[1] + yVecLocal[2]*i1_loc[2])

    return p, i1
end