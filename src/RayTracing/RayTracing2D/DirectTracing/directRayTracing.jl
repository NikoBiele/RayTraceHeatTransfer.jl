function direct_ray_tracing!(rtm::RayTracingDomain2D, rays_tot::P, nudge::G) where {G, P<:Integer}
    
    if rtm.spectral_mode != :grey
        println("Running direct ray tracing for $(rtm.n_spectral_bins) spectral bins")
        # Run direct ray tracing for each spectral bin
        for bin in 1:rtm.n_spectral_bins
            println("Processing spectral bin $bin/$(rtm.n_spectral_bins)")
            direct_ray_tracing_single_bin!(rtm, rays_tot, nudge, bin)
        end
    else
        println("Running direct ray tracing for grey extinction")
        direct_ray_tracing_single_bin!(rtm, rays_tot, nudge, 1)  # bin=1 for grey
    end

    update_scalar_temperatures_and_heat_sources_direct!(rtm)

end

function direct_ray_tracing_single_bin!(rtm::RayTracingDomain2D, rays_tot::P, nudge::G,
                                        spectral_bin::P) where {G, P<:Integer}

    # Prepare emitters
    emitters, total_energy = prepare_emitters(rtm, nudge, spectral_bin) # pass nudge to get the type G
    if total_energy == 0.0
        @warn "No emitters found for spectral bin $spectral_bin, skipping ray tracing"
        return
    end

    # Initialize counters for this spectral bin
    absorbed_count = [zeros(Int, length(coarse_face.subFaces)) for coarse_face in rtm.coarse_mesh]
    gas_emitted_count = [zeros(Int, length(coarse_face.subFaces)) for coarse_face in rtm.coarse_mesh]
    scattered_count = [zeros(Int, length(coarse_face.subFaces)) for coarse_face in rtm.coarse_mesh]
    wall_emitted_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subFaces] for coarse_face in rtm.coarse_mesh]
    reflected_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subFaces] for coarse_face in rtm.coarse_mesh]
    wall_absorbed_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subFaces] for coarse_face in rtm.coarse_mesh]

    nthreads = Threads.nthreads()
    
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
    progress = Progress(rays_tot, 1, "  Bin $spectral_bin ray tracing: ")
    completed_work = Threads.Atomic{Int}(0)
    
    # Thread-safe locks for writing to counters
    counter_locks = [Threads.SpinLock() for _ in 1:length(rtm.coarse_mesh)]
    
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
                # Check prescribed temperature for this spectral bin
                temp_value = fine_face.T_in_w[emitter.wall_index]
                if temp_value >= 0.0
                    local_wall_emitted_count[emitter.coarse_index][emitter.fine_index][emitter.wall_index] += 1
                end
            else
                fine_face = rtm.coarse_mesh[emitter.coarse_index].subFaces[emitter.fine_index]
                origin, direction = emit_volume_ray(fine_face, nudge)
                # Check prescribed temperature for this spectral bin
                temp_value = fine_face.T_in_g
                if temp_value >= 0.0
                    local_gas_emitted_count[emitter.coarse_index][emitter.fine_index] += 1
                end
            end
            
            # Use spectral ray tracing interface
            result = trace_single_ray(rtm, origin, direction, nudge, emitter.coarse_index, spectral_bin, 100000)
            
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
    
    update_spectral_results!(rtm, absorbed_count, gas_emitted_count, wall_emitted_count, 
                    reflected_count, scattered_count, wall_absorbed_count, 
                    total_energy, rays_tot, spectral_bin)
end