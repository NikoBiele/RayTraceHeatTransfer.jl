function emit_surface_ray(face::PolyFace2D{G}, wall_index::P, nudge::G) where {G, P<:Integer}
    # println("Emitting surface ray, type of G is $(typeof(G))")

    p1, p2 = face.vertices[wall_index], face.vertices[mod1(wall_index+1, length(face.vertices))]
    R = G(rand()) # Convert random number to the appropriate type G

    # Unit vectors in global coordinate system.
    xVecGlobal = Point2{G}(1.0, 0.0)
    yVecGlobal = Point2{G}(0.0, 1.0)

    # sample position to find emission point
    p = p1 + (p2 - p1) * R
    # nudge the point a tiny bit towards the midpoint, to ensure we are inside cell
    p = p + (face.midPoint - p) * nudge
    p = Point2{G}(p[1], p[2])

    # sample direction of ray (local coordinate system)
    i1_loc = lambertSample3D() # Convert the result to type G if needed
    
    # unit vectors in local coordinate system
    xVecLocal = normalize(p2-p1)
    yVecLocal = Point2{G}(-xVecLocal[2], xVecLocal[1])
    
    # rotate according to rotation matrix
    RotationMatrix = [dot(xVecGlobal, xVecLocal) dot(xVecGlobal, yVecLocal);
                dot(yVecGlobal, xVecLocal) dot(yVecGlobal, yVecLocal)]
    i1_init = RotationMatrix*i1_loc # emission direction (global coordinate system)
    i1 = Point2{G}(i1_init[1], i1_init[2])

    return p, i1
end