function directRayTracing!(rtm::RayTracingDomain2D, rays_tot::P, nudge::G, verbose::Bool,
                            seeds::Union{UnitRange{S},Vector{S}}, rngs::Vector{<:AbstractRNG}, nthreads::K) where {G,S<:Integer,P<:Integer,K<:Integer}
    
    if rtm.spectral_mode != :grey
        verbose && println("Running direct ray tracing for $(rtm.n_spectral_bins) spectral bins")
        # Run direct ray tracing for each spectral bin
        for bin in 1:rtm.n_spectral_bins
            verbose && println("Processing spectral bin $bin/$(rtm.n_spectral_bins)")
            directRayTracingSingleBin!(rtm, rays_tot, nudge, bin, seeds, rngs, nthreads; verbose=verbose)
        end
    else
        verbose && println("Running direct ray tracing for grey extinction")
        directRayTracingSingleBin!(rtm, rays_tot, nudge, 1, seeds, rngs, nthreads; verbose=verbose)  # bin=1 for grey
    end

    writeTemperaturesHeatSourcesDirect!(rtm; verbose=verbose)

end

function directRayTracingSingleBin!(rtm::RayTracingDomain2D, rays_tot::S, nudge::G,
                                        spectral_bin::P,
                                        seeds::Union{UnitRange{K},Vector{K}}, rngs::Vector{<:AbstractRNG},
                                        nthreads::P; verbose::Bool=true) where {G,S<:Integer,P<:Integer,K<:Integer}

    # Prepare emitters
    emitters, total_energy = prepareEmitters(rtm, nudge, spectral_bin) # pass nudge to get the type G
    if total_energy == 0.0
        @warn "No emitters found for spectral bin $spectral_bin, skipping ray tracing"
        return
    end

    # cumulative emitter energies, built once: an emitter is then picked with one random number and a binary search
    cdf = cumsum([e.energy for e in emitters])

    # Initialize counters for this spectral bin
    absorbed_count = [zeros(Int, length(coarse_face.subVolumes)) for coarse_face in rtm.coarse_mesh]
    gas_emitted_count = [zeros(Int, length(coarse_face.subVolumes)) for coarse_face in rtm.coarse_mesh]
    scattered_count = [zeros(Int, length(coarse_face.subVolumes)) for coarse_face in rtm.coarse_mesh]
    wall_emitted_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subVolumes] for coarse_face in rtm.coarse_mesh]
    reflected_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subVolumes] for coarse_face in rtm.coarse_mesh]
    wall_absorbed_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subVolumes] for coarse_face in rtm.coarse_mesh]
    
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
    
    # Progress tracking (only when verbose; assigned once so the threaded loop sees a concrete type)
    progress = verbose ? Progress(rays_tot; dt = 1, desc = "  Bin $spectral_bin ray tracing: ") : nothing
    completed_work = Threads.Atomic{Int}(0)
    
    # Thread-safe locks for writing to counters
    counter_locks = [Threads.SpinLock() for _ in 1:length(rtm.coarse_mesh)]
    
    @threads for tid in 1:nthreads
        ray_range = thread_assignments[tid]
        local_rng = rngs[tid] # Use thread ID as seed
        Random.seed!(rngs[tid], seeds[tid])
        path = Tuple{Int,Int,Int,Symbol}[]      # interactions of the current ray, reused for every ray
        
        # Local counters for this thread
        local_absorbed_count = [zeros(Int, length(coarse_face.subVolumes)) for coarse_face in rtm.coarse_mesh]
        local_gas_emitted_count = [zeros(Int, length(coarse_face.subVolumes)) for coarse_face in rtm.coarse_mesh]
        local_scattered_count = [zeros(Int, length(coarse_face.subVolumes)) for coarse_face in rtm.coarse_mesh]
        local_wall_emitted_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subVolumes] for coarse_face in rtm.coarse_mesh]
        local_reflected_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subVolumes] for coarse_face in rtm.coarse_mesh]
        local_wall_absorbed_count = [[zeros(Int, length(face.solidWalls)) for face in coarse_face.subVolumes] for coarse_face in rtm.coarse_mesh]

        for ray in ray_range
            # pick an emitter in proportion to its energy
            u = rand(local_rng, G) * cdf[end]
            emitter = emitters[min(searchsortedfirst(cdf, u), length(emitters))]
            
            if emitter.type == :surface
                fine_face = rtm.coarse_mesh[emitter.coarse_index].subVolumes[emitter.fine_index]
                origin, direction = emitSurfaceRay2D(fine_face, emitter.wall_index, nudge, local_rng)
                # Check prescribed temperature for this spectral bin
                temp_value = fine_face.T_in_w[emitter.wall_index]
                if temp_value >= 0.0
                    local_wall_emitted_count[emitter.coarse_index][emitter.fine_index][emitter.wall_index] += 1
                end
            else
                fine_face = rtm.coarse_mesh[emitter.coarse_index].subVolumes[emitter.fine_index]
                origin, direction = emitVolumeRay2D(fine_face, nudge, local_rng)
                # Check prescribed temperature for this spectral bin
                temp_value = fine_face.T_in_g
                if temp_value >= 0.0
                    local_gas_emitted_count[emitter.coarse_index][emitter.fine_index] += 1
                end
            end
            
            # Use spectral ray tracing interface
            result = traceSingleRay!(path, rtm, origin, direction, nudge, emitter.coarse_index, spectral_bin, 100_000, local_rng)
            
            if result !== nothing
                abs_coarse_index, abs_fine_index, abs_wall_index = result
                
                if abs_wall_index != 0                  # absorbed by a wall
                    local_wall_absorbed_count[abs_coarse_index][abs_fine_index][abs_wall_index] += 1
                else                                    # absorbed by the gas of that cell
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
            
            # Update progress in batches: one atomic update per 1000 rays instead of one per ray
            if verbose && (ray - first(ray_range) + 1) % 1000 == 0
                Threads.atomic_add!(completed_work, 1000)
                tid == 1 && update!(progress, completed_work[])
            end
        end
        
        # Merge local counters into global counters (thread-safe)
        for coarse_idx in 1:length(rtm.coarse_mesh)
            Threads.lock(counter_locks[coarse_idx]) do
                for fine_idx in 1:length(rtm.coarse_mesh[coarse_idx].subVolumes)
                    absorbed_count[coarse_idx][fine_idx] += local_absorbed_count[coarse_idx][fine_idx]
                    gas_emitted_count[coarse_idx][fine_idx] += local_gas_emitted_count[coarse_idx][fine_idx]
                    scattered_count[coarse_idx][fine_idx] += local_scattered_count[coarse_idx][fine_idx]
                    
                    for wall_idx in 1:length(rtm.coarse_mesh[coarse_idx].subVolumes[fine_idx].solidWalls)
                        wall_emitted_count[coarse_idx][fine_idx][wall_idx] += local_wall_emitted_count[coarse_idx][fine_idx][wall_idx]
                        reflected_count[coarse_idx][fine_idx][wall_idx] += local_reflected_count[coarse_idx][fine_idx][wall_idx]
                        wall_absorbed_count[coarse_idx][fine_idx][wall_idx] += local_wall_absorbed_count[coarse_idx][fine_idx][wall_idx]
                    end
                end
            end
        end
    end
    
    verbose && finish!(progress)
    
    updateSpectralResults!(rtm, absorbed_count, gas_emitted_count, wall_emitted_count, 
                    reflected_count, scattered_count, wall_absorbed_count, 
                    total_energy, rays_tot; spectral_bin=spectral_bin, verbose=verbose)
end