# Updated distToSurface function
function distToSurface(point::Point2{G}, direction::Point2{G}, face::PolyFace2D{G}) where G
    u = zeros(G, length(face.vertices))
    @simd for i in eachindex(face.vertices)
        j = mod1(i + 1, length(face.vertices))
        edge = face.vertices[j] - face.vertices[i]
        normal = -face.outwardNormals[i] # inward normal
        
        denominator = dot(direction, normal)
        if abs(denominator) < 1e-10  # Check for near-parallel rays
            u[i] = Inf
        else
            u[i] = dot(face.vertices[i] - point, normal) / denominator
        end
    end
    u[u .<= 0] .= Inf
    return findmin(u)
end