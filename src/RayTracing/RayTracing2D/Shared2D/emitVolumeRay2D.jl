function emitVolumeRay2D(face::PolyVolume2D{G}, nudge::G, rng::AbstractRNG) where {G}

    vertices = face.vertices
    if length(vertices) == 4
        A, B, C, D = vertices
        R_1, R_2 = rand(rng, G), rand(rng, G) # Convert to type G
        if rand(rng, G) < 0.5*(A[1]*(B[2]-C[2])+B[1]*(C[2]-A[2])+C[1]*(A[2]-B[2]))/face.volume
            sqrt_R1 = sqrt(R_1)
            point = (1-sqrt_R1)*A + sqrt_R1*(1-R_2)*B + sqrt_R1*R_2*C
        else
            sqrt_R1 = sqrt(R_1)
            point = (1-sqrt_R1)*C + sqrt_R1*(1-R_2)*D + sqrt_R1*R_2*A
        end
    else
        A, B, C = vertices
        R_1, R_2 = rand(rng, G), rand(rng, G) # Convert to type G
        sqrt_R1 = sqrt(R_1)
        point = (1-sqrt_R1)*A + sqrt_R1*(1-R_2)*B + sqrt_R1*R_2*C
    end

    # nudge the point a tiny bit towards the midpoint, to ensure we are inside cell
    point = point + (face.midPoint - point) * nudge
    point = Point2{G}(point[1], point[2])

    # uniform emission distribution
    theta = acos(1 - 2*rand(rng, G)) # Convert to type G
    phi = 2 * G(π) * rand(rng, G) # Convert to type G
    dir = Point2{G}(
        sin(theta)*cos(phi),
        cos(theta)
    )

    return point, dir
end