function find_nearest_face(faces::Vector{PolyFace2D{G,T}}, point::Point2{G}) where {G<:AbstractFloat, T}
    nearest_index = 1
    min_distance = Inf

    for (index, face) in enumerate(faces)
        distance = norm(point - face.midPoint)
        if distance < min_distance
            min_distance = distance
            nearest_index = index
        end
    end

    return nearest_index
end