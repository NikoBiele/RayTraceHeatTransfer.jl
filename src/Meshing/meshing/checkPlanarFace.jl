function checkPlanarFace(v::Vector{Point3{G}}, i::Int, tol::G) where {G}
    length(v) < 4 && return
    n = normalize(cross(v[2] - v[1], v[3] - v[1]))
    scale = maximum(norm(v[k] - v[1]) for k in 2:4)
    d = abs(dot(v[4] - v[1], n))
    if d > tol * scale
        error("Face $i is non-planar (out-of-plane deviation $(d/scale) of face size). " *
              "The mesher projects onto the plane of the first three vertices, so the " *
              "fourth vertex would be silently moved. Split this face into two triangles.")
    end
end