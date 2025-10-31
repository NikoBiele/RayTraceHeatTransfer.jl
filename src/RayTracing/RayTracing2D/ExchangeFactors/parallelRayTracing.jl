function parallel_ray_tracing_optimized(rtm::RayTracingDomain2D, rays_total::P, 
                                        nudge::G) where {P<:Integer, G}

    surface_mapping, volume_mapping, num_surfaces, num_volumes = create_index_mapping(rtm, rays_total)
    num_emitters = num_surfaces + num_volumes
    rays_per_emitter = div(rays_total, num_emitters)

    # Determine spectral configuration
    n_bins = rtm.n_spectral_bins
    is_spectral = rtm.spectral_mode != :grey
    
    if rtm.spectral_mode == :spectral_variable
        # Variable spectral - need separate F matrix for each bin
        println("Computing $n_bins separate F matrices for variable spectral extinction")
        F_raw_vector = Matrix{G}[]
        # F_raw_uncertain_vector = Matrix{Measurements.Measurement{G}}[]
        
        for bin in 1:n_bins
            println("Computing F matrix for spectral bin $bin/$n_bins")
            F_raw_bin = #, F_raw_uncertain_bin = 
                compute_F_matrix_for_bin(
                rtm, rays_per_emitter, nudge, bin,
                surface_mapping, volume_mapping, num_surfaces, num_volumes, num_emitters
            )
            push!(F_raw_vector, F_raw_bin)
            # push!(F_raw_uncertain_vector, F_raw_uncertain_bin)
        end
        
        return F_raw_vector, rays_per_emitter # F_raw_uncertain_vector,
        
    else
        # Grey or uniform spectral - single F matrix works for all bins
        if is_spectral
            println("Computing single F matrix for uniform spectral extinction ($n_bins bins)")
        else
            println("Computing single F matrix for grey extinction")
        end
        
        F_raw = compute_F_matrix_for_bin(
            rtm, rays_per_emitter, nudge, 1,  # Use bin 1 (doesn't matter for uniform)
            surface_mapping, volume_mapping, num_surfaces, num_volumes, num_emitters
        )
        
        return F_raw, rays_per_emitter # F_raw_uncertain,
    end
end

function compute_F_matrix_for_bin(rtm::RayTracingDomain2D, rays_per_emitter::P, 
                                 nudge::G, spectral_bin::P,
                                 surface_mapping, volume_mapping, num_surfaces, num_volumes, num_emitters) where {P<:Integer, G}
    
    # Pre-allocate result matrix
    F_counts = zeros(Int, num_emitters, num_emitters)
    
    # Create combined emitter list
    all_emitters = Vector{Tuple{Any, Int}}()
    for emitter_key in keys(surface_mapping)
        global_idx = surface_mapping[emitter_key]
        push!(all_emitters, (emitter_key, global_idx))
    end
    for emitter_key in keys(volume_mapping)
        global_idx = num_surfaces + volume_mapping[emitter_key]
        push!(all_emitters, (emitter_key, global_idx))
    end
    sort!(all_emitters, by=x->x[2])
    
    nthreads = Threads.nthreads()
    println("  Using $nthreads threads for spectral bin $spectral_bin")
    
    # Divide emitters among threads
    emitters_per_thread = div(num_emitters, nthreads)
    remainder = num_emitters % nthreads
    
    thread_assignments = Vector{UnitRange{Int}}(undef, nthreads)
    start_idx = 1
    for tid in 1:nthreads
        thread_size = emitters_per_thread + (tid <= remainder ? 1 : 0)
        end_idx = start_idx + thread_size - 1
        thread_assignments[tid] = start_idx:end_idx
        start_idx = end_idx + 1
    end
    
    # Progress tracking
    progress = Progress(num_emitters, 1, "  Bin $spectral_bin progress: ")
    completed_work = Threads.Atomic{Int}(0)
    
    # Thread-safe locks for writing to F_counts
    row_locks = [Threads.SpinLock() for _ in 1:num_emitters]
    
    # Parallel ray tracing for this spectral bin
    @threads for tid in 1:nthreads
        emitter_range = thread_assignments[tid]
        
        for global_emitter_idx in emitter_range
            emitter_key, global_idx = all_emitters[global_emitter_idx]
            
            # Determine emitter type and process rays
            is_surface_emitter = global_idx <= num_surfaces
            
            # Temporary storage for this emitter's results
            temp_row = zeros(Int, num_emitters)
            
            if is_surface_emitter
                coarse_index, fine_index, wall_index = emitter_key
                face = rtm.fine_mesh[coarse_index][fine_index]
                
                for _ in 1:rays_per_emitter
                    p_emit, dir_emit = emit_surface_ray(face, wall_index, nudge)
                    # Pass spectral bin to trace_ray
                    result = trace_ray(rtm, p_emit, dir_emit, nudge, coarse_index, spectral_bin)
                    
                    if result !== nothing
                        absorber_index = get_global_index(surface_mapping, volume_mapping, num_surfaces, result...)
                        if absorber_index != -1
                            temp_row[absorber_index] += 1
                        end
                    end
                end
            else
                coarse_index, fine_index = emitter_key
                face = rtm.fine_mesh[coarse_index][fine_index]
                
                for _ in 1:rays_per_emitter
                    p_emit, dir_emit = emit_volume_ray(face, nudge)
                    # Pass spectral bin to trace_ray
                    result = trace_ray(rtm, p_emit, dir_emit, nudge, coarse_index, spectral_bin)
                    
                    if result !== nothing
                        absorber_index = get_global_index(surface_mapping, volume_mapping, num_surfaces, result...)
                        if absorber_index != -1
                            temp_row[absorber_index] += 1
                        end
                    end
                end
            end
            
            # Write results to final matrix (thread-safe)
            Threads.lock(row_locks[global_idx]) do
                for j in 1:num_emitters
                    F_counts[global_idx, j] += temp_row[j]
                end
            end
            
            # Update progress
            Threads.atomic_add!(completed_work, 1)
            if tid == 1
                update!(progress, completed_work[])
            end
        end
    end
    finish!(progress)
    
    # Convert to exchange factors
    F_raw = zeros(G, num_emitters, num_emitters)
    # F_raw_uncertain = zeros(G, num_emitters, num_emitters)

    progress_f = Progress(num_emitters^2, 1, "  Calculating F for bin $spectral_bin: ")
    calculated = 0
    
    for i in 1:num_emitters
        for j in 1:num_emitters
            if F_counts[i, j] > 0
                if G <: Measurement
                    F_raw[i, j] = (F_counts[i, j] / rays_per_emitter) ± (sqrt(F_counts[i, j]) / rays_per_emitter)
                else
                    F_raw[i, j] = F_counts[i, j] / rays_per_emitter
                end
            end
            calculated += 1
            if calculated % 1000 == 0
                update!(progress_f, calculated)
            end
        end
    end
    finish!(progress_f)

    return F_raw
end