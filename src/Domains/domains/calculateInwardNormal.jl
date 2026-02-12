function calculateInwardNormal(p1::Point2{G}, p2::Point2{G}, midpoint::Point2{G}) where G
    edge = p2 - p1
    normal = normalize(Point2{G}(edge[2], -edge[1])) # normal
    
    # Check if normal is pointing outward
    wall_midpoint = (p1 + p2) / 2
    if dot(normal, wall_midpoint - midpoint) < 0
        normal = -normal  # Flip if pointing outward
    end
    
    return normal
end

function calculateInwardNormal(p1::Point3{G}, p2::Point3{G}, p3::Point3{G}, midpoint::Point3{G}) where {G}
    edge1 = p2 - p1
    edge2 = p3 - p1
    normal = normalize(cross(edge1, edge2))
    
    face_midpoint = (p1 + p2 + p3) / 3
    if dot(normal, midpoint - face_midpoint) < 0
        normal = -normal
    end
    
    return normal
end