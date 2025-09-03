function exchangeFactors3D!(superFaces::Vector{Face3D})

    # calculate the view factors
    num_faces = length(superFaces)
    num_subFaces = length(superFaces[1].subFaces)
    elements = num_faces * num_subFaces
    F = zeros(elements, elements)
    areas_linear = Vector{Float64}(undef, elements)

    # emitter count
    emitter_count = 0
    for p in 1:num_faces # for each coarse face
        for j in 1:num_subFaces

            # increment the emitter count
            emitter_count += 1

            # set the emitting face
            if length(superFaces[p].subFaces[j].vertices) == 3
                face1 = Matrix{Float64}(undef, 3, 3)
            else
                face1 = Matrix{Float64}(undef, 4, 3)
            end
            for i = 1:size(face1,1)
                face1[i, :] = superFaces[p].subFaces[j].vertices[i]
            end

            # absorber count
            absorber_count = 0
            # find the absorbing face
            for k in 1:num_faces
                for l in 1:num_subFaces

                    # increment the absorber count
                    absorber_count += 1

                    # set the absorbing face
                    if length(superFaces[k].subFaces[l].vertices) == 3
                        face2 = Matrix{Float64}(undef, 3, 3)
                    else
                        face2 = Matrix{Float64}(undef, 4, 3)
                    end
                    for i = 1:size(face2,1)    
                        face2[i, :] = superFaces[k].subFaces[l].vertices[i]
                    end

                    # viewFactor
                    result = viewFactor(face1, face2)
                    if emitter_count != absorber_count # should avoid the 'emitter == absorber' case
                        F[emitter_count, absorber_count] = occursin("NaN", string(result[1])) ? 0.0 : result[1]
                    else
                        F[emitter_count, absorber_count] = 0.0
                    end
                    superFaces[p].subFaces[j].area = result[3]
                    areas_linear[emitter_count] = result[3]
                end
            end
        end
    end

    # F, _ = smooth_exchange_factors_ultimate!(F, areas_linear, nothing, GasProperties(1.0, 1.0), 1, false; max_iterations=200, tolerance=eps(Float64))

    # For 3D case with only surface elements (no participating media)
    surface_areas = areas_linear
    volumes = Float64[]  # Empty vector since no volume elements
    volume_betas = Float64[]  # Empty vector since no volume elements

    F, _ = smooth_exchange_factors_ultimate!(F, surface_areas, volumes, volume_betas, 
                                        1, false; max_iterations=200, tolerance=eps(Float64))

    return F
    
end