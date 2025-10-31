function exchangeFactors3D!(superFaces::Vector{PolyFace3D})
    # Calculate dimensions
    num_faces = length(superFaces)
    num_subFaces = length(superFaces[1].subFaces)
    elements = num_faces * num_subFaces
    
    # Pre-allocate results
    F = zeros(elements, elements)
    areas_linear = Vector{Float64}(undef, elements)
    
    # Convert 2D indices to linear index
    @inline function to_linear(face_idx, subface_idx)
        return (face_idx - 1) * num_subFaces + subface_idx
    end
    
    # Convert linear index back to 2D indices  
    @inline function from_linear(linear_idx)
        face_idx = div(linear_idx - 1, num_subFaces) + 1
        subface_idx = mod(linear_idx - 1, num_subFaces) + 1
        return face_idx, subface_idx
    end
    
    # Parallelize over emitters (outer loop)
    @threads for emitter_linear in 1:elements
        emitter_face, emitter_sub = from_linear(emitter_linear)
        
        # Set up emitting face geometry
        vertices1 = superFaces[emitter_face].subFaces[emitter_sub].vertices
        face1 = Matrix{Float64}(undef, length(vertices1), 3)
        for i in 1:length(vertices1)
            face1[i, :] = vertices1[i]
        end
        
        # Loop over all absorbers for this emitter
        for absorber_linear in 1:elements
            if emitter_linear == absorber_linear
                F[emitter_linear, absorber_linear] = 0.0
                continue
            end
            
            absorber_face, absorber_sub = from_linear(absorber_linear)
            
            # Set up absorbing face geometry
            vertices2 = superFaces[absorber_face].subFaces[absorber_sub].vertices
            face2 = Matrix{Float64}(undef, length(vertices2), 3)
            for i in 1:length(vertices2)
                face2[i, :] = vertices2[i]
            end
            
            # Calculate view factor
            result = viewFactor(face1, face2)
            F[emitter_linear, absorber_linear] = occursin("NaN", string(result[1])) ? 0.0 : result[1]
        end
        
        # Store area (only needs to be done once per emitter)
        result = viewFactor(face1, face1)  # Get area from self-calculation
        areas_linear[emitter_linear] = result[3]
        superFaces[emitter_face].subFaces[emitter_sub].area = result[3]
    end
    
    # Apply smoothing
    surface_areas = areas_linear
    volumes = Float64[]
    volume_betas = Float64[]
    
    F = smooth_exchange_factors_ultimate!(F, surface_areas, volumes, volume_betas, 
                                           1; max_iterations=200, tolerance=eps(Float64))
    
    return F
end