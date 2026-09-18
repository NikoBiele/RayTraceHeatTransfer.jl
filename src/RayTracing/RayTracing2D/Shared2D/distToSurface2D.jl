function distToSurface2D(point::Point2{G}, direction::Point2{G}, face::PolyVolume2D{G}) where G
    verts = face.vertices
    normals = face.inwardNormals
    n = length(verts)
    u_min = G(Inf)
    i_min = 0
    @inbounds for i in 1:n
        normal = normals[i]
        denominator = dot(direction, normal)
        abs(denominator) < 1e-10 && continue                 # parallel to this edge
        u = dot(verts[i] - point, normal) / denominator
        if u > 0 && u < u_min
            u_min = u
            i_min = i
        end
    end
    return u_min, i_min
end