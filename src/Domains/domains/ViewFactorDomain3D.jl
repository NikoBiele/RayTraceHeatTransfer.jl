


# Updated constructor with spectral support
function ViewFactorDomain3D(points::Matrix{G}, faces::Matrix{P}, Ndims::P,
                       q_in_w::Vector{G}, T_in_w::Vector{G}, 
                       epsilon::Union{Vector{G}, Vector{Vector{G}}}) where {G, P<:Integer}
    
    # Determine spectral mode from epsilon
    is_spectral = isa(epsilon[1], Vector)
    if is_spectral
        if std([std(epsilon[i]) for i in 1:size(epsilon, 1)]) > 1e-6
            spectral_mode = :spectral_variable
        else
            spectral_mode = :spectral_uniform
        end
    else
        spectral_mode = :grey
    end
    n_bins = is_spectral ? length(epsilon[1]) : 1
    
    # Calculate domain midpoint
    domain_mid = sum(points, dims=1)/size(points, 1)    
    domain_midpoint = Point3{G}(domain_mid[1], domain_mid[2], domain_mid[3])
    
    # Create super faces
    superFaces = PolyFace3D[]
    for (i, face_rows) in enumerate(eachrow(faces))
        if length(face_rows) == 4
            points3d = [Point3{G}(points[face_rows[j],:]) for j in 1:4]
        else
            points3d = [Point3{G}(points[face_rows[j],:]) for j in 1:3]
        end
        
        push!(superFaces, PolyFace3D(points3d, true, domain_midpoint, epsilon[i], q_in_w[i], T_in_w[i]))
    end
    
    # Mesh the faces
    mesh3D = meshFaces(points, faces, Ndims)
    num_faces = length(mesh3D)
    num_points = length(mesh3D[1][1])
    
    first_epsilon = nothing
    uniform_epsilon = true
    for i in 1:num_faces
        superFaces[i].subFaces = PolyFace3D[]
        
        # First pass: create all subfaces with q=0
        for j in 1:num_points
            p1 = Point3{G}(mesh3D[i][1][j]...)
            p2 = Point3{G}(mesh3D[i][2][j]...)
            p3 = Point3{G}(mesh3D[i][3][j]...)
            p4 = Point3{G}(mesh3D[i][4][j]...)
            
            if isapprox(p3, p4, atol=1e-5) 
                push!(superFaces[i].subFaces, 
                    PolyFace3D([p1, p2, p3], true, domain_midpoint, epsilon[i], 0.0, T_in_w[i]))
            else
                push!(superFaces[i].subFaces, 
                    PolyFace3D([p1, p2, p3, p4], true, domain_midpoint, epsilon[i], 0.0, T_in_w[i]))
            end
        end
        
        # Second pass: distribute flux proportional to area
        total_area = sum(subface.area for subface in superFaces[i].subFaces)
        
        for subface in superFaces[i].subFaces
            # Each subface gets flux proportional to its area fraction
            subface.q_in_w = q_in_w[i] * (subface.area / total_area)
            if first_epsilon === nothing && isa(subface.epsilon, Vector)
                first_epsilon = subface.epsilon[1]
            elseif first_epsilon === nothing && !isa(subface.epsilon, Vector)
                first_epsilon = subface.epsilon
            else
                for bin in 1:n_bins
                    if !isapprox(subface.epsilon[bin], first_epsilon, atol=1e-5)
                        uniform_epsilon = false
                    end
                end
            end
        end
    end

    # Return with superFaces as facesMesh (this is the variable name expected)
    number_of_elements = size(faces)[1]
    F_raw = Matrix{G}(undef, number_of_elements, number_of_elements)
    F_smooth = Matrix{G}(undef, number_of_elements, number_of_elements)
    return ViewFactorDomain3D{G, P}(points, faces, Ndims, superFaces, F_raw, F_smooth, 
                                spectral_mode, n_bins, nothing, nothing, uniform_epsilon)
end
    
function (vfd::ViewFactorDomain3D)(; parallel::Bool=true, tol::G=1e-15) where G

    # Calculate view factors (wavelength-independent!)
    println("Computing view factors (geometry only, wavelength-independent)...")
    F_raw, F_smooth = enclosureViewFactors3D(vfd.facesMesh, parallel, tol)
    
    vfd.F_raw = F_raw
    vfd.F_smooth = F_smooth

    return nothing
end