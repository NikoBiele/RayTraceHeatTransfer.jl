function parallel_ray_tracing_optimized(rtm::RayTracingMeshOptim, rays_total::P, 
                                         uncertain::Bool, nudge=Float64(eps(Float32))) where {P<:Integer}

    surface_mapping, volume_mapping, num_surfaces, num_volumes = create_index_mapping(rtm, rays_total)
    num_emitters = num_surfaces + num_volumes
    rays_per_emitter = div(rays_total, num_emitters)

    # Pre-allocate final result matrix
    F_counts = zeros(P, num_emitters, num_emitters)
    
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
    println("Using $nthreads threads for in-place ray tracing")
    
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
    
    println("Memory usage: $(sizeof(P) * num_emitters^2 / 1e9) GB (single matrix)")
    
    # Progress tracking
    progress = Progress(num_emitters, 1, "Ray tracing progress: ")
    completed_work = Threads.Atomic{Int}(0)
    
    # Thread-safe locks for writing to F_counts
    row_locks = [Threads.SpinLock() for _ in 1:num_emitters]
    
    # Parallel ray tracing - write directly to final matrix
    println("Starting in-place ray tracing...")
    @threads for tid in 1:nthreads
        emitter_range = thread_assignments[tid]
        
        for global_emitter_idx in emitter_range
            emitter_key, global_idx = all_emitters[global_emitter_idx]
            
            # Determine emitter type and process rays
            is_surface_emitter = global_idx <= num_surfaces
            
            # Temporary storage for this emitter's results
            temp_row = zeros(P, num_emitters)
            
            if is_surface_emitter
                coarse_index, fine_index, wall_index = emitter_key
                face = rtm.fine_mesh[coarse_index][fine_index]
                
                for _ in 1:rays_per_emitter
                    p_emit, dir_emit = emit_surface_ray(face, wall_index, nudge)
                    # Updated: removed extinction parameters
                    result = trace_ray(rtm, p_emit, dir_emit, nudge, coarse_index)
                    
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
                    # Updated: removed extinction parameters
                    result = trace_ray(rtm, p_emit, dir_emit, nudge, coarse_index)
                    
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
    println("In-place ray tracing complete.")
    
    # Convert to exchange factors
    # Determine the appropriate type for F_raw
    if P <: Measurement
        NumericType = P.parameters[1]  # Get underlying numeric type from Measurement
    else
        NumericType = Float64  # Default to Float64 for integer types
    end
    
    F_raw = zeros(NumericType, num_emitters, num_emitters)
    
    if uncertain
        F_raw_uncertain = zeros(Measurement{NumericType}, num_emitters, num_emitters)
    else
        F_raw_uncertain = zeros(NumericType, num_emitters, num_emitters)
    end

    progress_f = Progress(num_emitters^2, 1, "Calculating F, progress: ")
    calculated = 0
    
    for i in 1:num_emitters
        for j in 1:num_emitters
            if F_counts[i, j] > zero(P)
                F_raw[i, j] = F_counts[i, j] / rays_per_emitter
            end
            if uncertain && F_counts[i, j] > zero(P)
                F_raw_uncertain[i, j] = (F_counts[i, j] / rays_per_emitter) ± (sqrt(F_counts[i, j]) / rays_per_emitter)
            end
            calculated += 1
            if calculated % 1000 == 0
                update!(progress_f, calculated)
            end
        end
    end
    finish!(progress_f)

    return F_raw, F_raw_uncertain, rays_per_emitter
end

# In-place accumulation (zero additional memory)
function sum_thread_results_inplace!(thread_results::Vector{Matrix{T}}) where T
    nthreads = length(thread_results)
    if nthreads == 1
        return thread_results[1]
    end
    
    n = size(thread_results[1], 1)
    
    # Accumulate everything into the first matrix
    result = thread_results[1]
    
    @threads for i in 1:n
        @inbounds for thread_id in 2:nthreads
            for j in 1:n
                result[i, j] += thread_results[thread_id][i, j]
            end
        end
    end
    
    return result
end