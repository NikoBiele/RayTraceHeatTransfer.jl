# Check reciprocity
function check_reciprocity(F::Matrix{P}, mesh::RayTracingMeshOptim, surface_mapping::Dict, volume_mapping::Dict, gas::GasProperties, tol=1e-12) where {P}
    num_emitters = size(F, 1)
    reciprocity_satisfied = true
    max_relative_error = 0.0
    violating_pairs = Tuple{Int, Int}[]

    num_surfaces = length(surface_mapping)
    num_volumes = length(volume_mapping)

    # Create arrays for surface areas and volumes
    surface_areas = zeros(P, num_surfaces)
    volumes = zeros(P, num_volumes)

    # Populate surface areas and volumes
    for (coarse_index, coarse_face) in enumerate(mesh.coarse_mesh)
        for (fine_index, fine_face) in enumerate(mesh.fine_mesh[coarse_index])
            for (wall_index, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    emitter_index = get(surface_mapping, (coarse_index, fine_index, wall_index), 0)
                    if emitter_index > 0
                        surface_areas[emitter_index] = fine_face.area[wall_index]
                    end
                end
            end
            emitter_index = get(volume_mapping, (coarse_index, fine_index), 0)
            if emitter_index > 0
                volumes[emitter_index] = fine_face.volume
            end
        end
    end

    # Check reciprocity
    for i in 1:num_emitters
        for j in (i+1):num_emitters
            left = 0.0
            right = 0.0

            if i <= num_surfaces && j <= num_surfaces
                # Surface-to-surface
                left = F[i, j] * surface_areas[i]
                right = F[j, i] * surface_areas[j]
            elseif i <= num_surfaces && j > num_surfaces
                # Surface-to-volume
                left = F[i, j] * surface_areas[i]
                right = F[j, i] * 4 * gas.beta * volumes[j - num_surfaces]
            elseif i > num_surfaces && j <= num_surfaces
                # Volume-to-surface
                left = F[i, j] * 4 * gas.beta * volumes[i - num_surfaces]
                right = F[j, i] * surface_areas[j]
            elseif i > num_surfaces && j > num_surfaces
                # Volume-to-volume
                left = 4 * gas.beta * F[i, j] * volumes[i - num_surfaces]
                right = 4 * gas.beta * F[j, i] * volumes[j - num_surfaces]
            end

            if !isapprox(left, right, atol=tol)
                reciprocity_satisfied = false
                relative_error = abs(left - right) / max(abs(left), abs(right))
                max_relative_error = max(max_relative_error, relative_error)
                push!(violating_pairs, (i, j))
                
                println("Debug: Violation at (i=$i, j=$j), left=$left, right=$right, relative_error=$relative_error")
                
                # Optionally, limit the number of violating pairs reported
                if length(violating_pairs) >= 10
                    break
                end
            end
        end
        if !reciprocity_satisfied && length(violating_pairs) >= 10
            break
        end
    end

    if reciprocity_satisfied
        println("Reciprocity is satisfied for all pairs within the specified tolerance.")
    else
        println("Reciprocity is not satisfied. Maximum relative error: $max_relative_error")
        println("First few violating pairs:")
        for (i, j) in violating_pairs[1:min(5, length(violating_pairs))]
            println("($i, $j)")
        end
    end

    return reciprocity_satisfied, max_relative_error, violating_pairs
end
