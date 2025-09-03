function direct_ray_tracing!(rtm::RayTracingMeshOptim, rays_tot::Int, nudge=Float64(eps(Float32)))

    emitters, total_energy = prepare_emitters(rtm)

    # gas interaction counters
    absorbed_count = [zeros(Int, length(coarse_face.subFaces)) for coarse_face in rtm.coarse_mesh]
    gas_emitted_count = [zeros(Int, length(coarse_face.subFaces)) for coarse_face in rtm.coarse_mesh]
    scattered_count = [zeros(Int, length(coarse_face.subFaces)) for coarse_face in rtm.coarse_mesh]
    # wall interaction counters
    wall_emitted_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subFaces] for coarse_face in rtm.coarse_mesh]
    reflected_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subFaces] for coarse_face in rtm.coarse_mesh]
    wall_absorbed_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subFaces] for coarse_face in rtm.coarse_mesh]

    nthreads = Threads.nthreads()
    println("Using $nthreads threads for direct ray tracing")
    
    # Divide rays among threads
    rays_per_thread = div(rays_tot, nthreads)
    remainder = rays_tot % nthreads
    
    thread_assignments = Vector{UnitRange{Int}}(undef, nthreads)
    start_idx = 1
    for tid in 1:nthreads
        thread_size = rays_per_thread + (tid <= remainder ? 1 : 0)
        end_idx = start_idx + thread_size - 1
        thread_assignments[tid] = start_idx:end_idx
        start_idx = end_idx + 1
    end
    
    # Progress tracking
    progress = Progress(rays_tot, 1, "Ray tracing progress: ")
    completed_work = Threads.Atomic{Int}(0)
    
    # Thread-safe locks for writing to counters
    counter_locks = [Threads.SpinLock() for _ in 1:length(rtm.coarse_mesh)]
    
    println("Starting direct ray tracing...")
    @threads for tid in 1:nthreads
        ray_range = thread_assignments[tid]
        local_rng = MersenneTwister(tid)  # Use thread ID as seed
        
        # Local counters for this thread
        local_absorbed_count = [zeros(Int, length(coarse_face.subFaces)) for coarse_face in rtm.coarse_mesh]
        local_gas_emitted_count = [zeros(Int, length(coarse_face.subFaces)) for coarse_face in rtm.coarse_mesh]
        local_scattered_count = [zeros(Int, length(coarse_face.subFaces)) for coarse_face in rtm.coarse_mesh]
        local_wall_emitted_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subFaces] for coarse_face in rtm.coarse_mesh]
        local_reflected_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subFaces] for coarse_face in rtm.coarse_mesh]
        local_wall_absorbed_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subFaces] for coarse_face in rtm.coarse_mesh]

        for ray in ray_range
            emitter = sample(local_rng, emitters, Weights([e.energy for e in emitters]))
            
            if emitter.type == :surface
                fine_face = rtm.coarse_mesh[emitter.coarse_index].subFaces[emitter.fine_index]
                origin, direction = emit_surface_ray(fine_face, emitter.wall_index, nudge)
                # Only count initial wall emissions for prescribed temperature walls
                if fine_face.T_in_w[emitter.wall_index] >= 0.0
                    local_wall_emitted_count[emitter.coarse_index][emitter.fine_index][emitter.wall_index] += 1
                end
            else
                fine_face = rtm.coarse_mesh[emitter.coarse_index].subFaces[emitter.fine_index]
                origin, direction = emit_volume_ray(fine_face, nudge)
                # Only count initial gas emissions for prescribed temperature gas
                if fine_face.T_in_g >= 0.0
                    local_gas_emitted_count[emitter.coarse_index][emitter.fine_index] += 1
                end
            end
            
            # Use new ray tracing interface
            result = trace_single_ray(rtm, origin, direction, nudge, emitter.coarse_index, 100000)
            
            if result !== nothing
                absorption_type, abs_coarse_index, abs_fine_index, abs_wall_index, path = result
                
                if absorption_type == :surface
                    local_wall_absorbed_count[abs_coarse_index][abs_fine_index][abs_wall_index] += 1
                elseif absorption_type == :gas
                    local_absorbed_count[abs_coarse_index][abs_fine_index] += 1
                end
                
                for (coarse_index, fine_index, wall_index, interaction_type) in path
                    if interaction_type == :reflection
                        local_reflected_count[coarse_index][fine_index][wall_index] += 1
                    elseif interaction_type == :scattering
                        local_scattered_count[coarse_index][fine_index] += 1
                    elseif interaction_type == :reemission
                        if wall_index == 0
                            # Gas reemission - count both absorption and emission
                            local_absorbed_count[coarse_index][fine_index] += 1
                            local_gas_emitted_count[coarse_index][fine_index] += 1
                        else
                            # Wall reemission - count both absorption and emission
                            local_wall_absorbed_count[coarse_index][fine_index][wall_index] += 1
                            local_wall_emitted_count[coarse_index][fine_index][wall_index] += 1
                        end
                    end
                end
            end
            
            # Update progress
            Threads.atomic_add!(completed_work, 1)
            if tid == 1 && (ray - ray_range.start + 1) % 100 == 0
                update!(progress, completed_work[])
            end
        end
        
        # Merge local counters into global counters (thread-safe)
        for coarse_idx in 1:length(rtm.coarse_mesh)
            Threads.lock(counter_locks[coarse_idx]) do
                for fine_idx in 1:length(rtm.coarse_mesh[coarse_idx].subFaces)
                    absorbed_count[coarse_idx][fine_idx] += local_absorbed_count[coarse_idx][fine_idx]
                    gas_emitted_count[coarse_idx][fine_idx] += local_gas_emitted_count[coarse_idx][fine_idx]
                    scattered_count[coarse_idx][fine_idx] += local_scattered_count[coarse_idx][fine_idx]
                    
                    for wall_idx in 1:length(rtm.coarse_mesh[coarse_idx].subFaces[fine_idx].solidWalls)
                        wall_emitted_count[coarse_idx][fine_idx][wall_idx] += local_wall_emitted_count[coarse_idx][fine_idx][wall_idx]
                        reflected_count[coarse_idx][fine_idx][wall_idx] += local_reflected_count[coarse_idx][fine_idx][wall_idx]
                        wall_absorbed_count[coarse_idx][fine_idx][wall_idx] += local_wall_absorbed_count[coarse_idx][fine_idx][wall_idx]
                    end
                end
            end
        end
    end
    
    finish!(progress)
    println("Direct ray tracing complete.")

    update_heat_source!(rtm, absorbed_count, gas_emitted_count, wall_emitted_count, 
                    reflected_count, scattered_count, wall_absorbed_count, 
                    total_energy, rays_tot)
end
