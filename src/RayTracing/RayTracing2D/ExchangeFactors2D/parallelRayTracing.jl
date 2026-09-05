function parallelRayTracing(rtm::RayTracingDomain2D, rays_total::S, 
                                        nudge::G, verbose::Bool,
                                        seeds::Union{UnitRange{P},Vector{P}}, rngs::Vector{<:AbstractRNG},
                                        nthreads::K;
                                        rec=nothing) where {S<:Integer, P<:Integer, K<:Integer, G}

    surface_mapping, volume_mapping, num_surfaces, num_volumes = createIndexMapping(rtm, rays_total)
    num_emitters = num_surfaces + num_volumes
    rays_per_emitter = div(rays_total, num_emitters)

    # Determine spectral configuration
    n_bins = rtm.n_spectral_bins
    is_spectral = rtm.spectral_mode != :grey
    
    if rtm.spectral_mode == :spectral_variable && !rtm.surfaces_only
        # Variable spectral - need separate F matrix for each bin
        verbose && println("Computing $n_bins separate F matrices for variable spectral extinction")
        F_raw_vector = Vector{AbstractMatrix}(undef, n_bins)
                
        groups, reps, nonuniform = group_uniform_bins(rtm.uniform_across_bin)

        # loop over nonuniform
        for bin in nonuniform
            verbose && println("Computing F matrix for nonuniform spectral bin $bin/$n_bins")
            F_raw_bin = computeExchangeFactorsBin(
                rtm, rays_per_emitter, nudge, bin,
                surface_mapping, volume_mapping, num_surfaces,
                num_volumes, num_emitters, verbose, rec, seeds, rngs, nthreads
            )
            F_raw_vector[bin] = F_raw_bin
        end

        # loop over groups
        for (k, idx_group) in pairs(groups)
            representative_bin = first(idx_group)
            verbose && println("Computing F matrix for uniform spectral bin $representative_bin/$n_bins")
            F_raw_bin = computeExchangeFactorsBin(
                rtm, rays_per_emitter, nudge, representative_bin,
                surface_mapping, volume_mapping, num_surfaces,
                num_volumes, num_emitters, verbose, rec, seeds, rngs, nthreads
            )
            for j in idx_group
                F_raw_vector[j] = F_raw_bin
            end
        end

        return F_raw_vector
        
    else
        # Grey or uniform spectral or surfaces only - single F matrix works for all bins
        if is_spectral
            verbose && println("Computing single F matrix for uniform spectral extinction ($n_bins bins)")
        else
            verbose && println("Computing single F matrix for grey extinction")
        end
        
        F_raw = computeExchangeFactorsBin(
            rtm, rays_per_emitter, nudge, 1,  # Use bin 1 (doesn't matter for uniform)
            surface_mapping, volume_mapping, num_surfaces,
            num_volumes, num_emitters, verbose, rec, seeds, rngs, nthreads
        )
        
        if rtm.surfaces_only
            return F_raw[1:length(rtm.surface_mapping),1:length(rtm.surface_mapping)]
        else
            return F_raw
        end
    end
end

function computeExchangeFactorsBin(rtm::RayTracingDomain2D, rays_per_emitter::S,
                                 nudge::G, spectral_bin::P,
                                 surface_mapping, volume_mapping, num_surfaces,
                                 num_volumes, num_emitters, verbose, rec,
                                 seeds::Union{UnitRange{K},Vector{K}}, rngs::Vector{<:AbstractRNG},
                                 nthreads::P) where {S<:Integer, P<:Integer, K<:Integer, G}

    # Combined, globally-sorted emitter list
    all_emitters = Vector{Tuple{Any, Int}}()
    for k in keys(surface_mapping)
        push!(all_emitters, (k, surface_mapping[k]))
    end
    for k in keys(volume_mapping)
        push!(all_emitters, (k, num_surfaces + volume_mapping[k]))
    end
    sort!(all_emitters, by = x -> x[2])

    verbose && println("  Using $nthreads threads for spectral bin $spectral_bin")

    # Partition emitters across threads — disjoint ranges, so each thread is independent
    emitters_per_thread = div(num_emitters, nthreads)
    remainder = num_emitters % nthreads
    thread_assignments = Vector{UnitRange{Int}}(undef, nthreads)
    start_idx = 1
    for tid in 1:nthreads
        thread_size = emitters_per_thread + (tid <= remainder ? 1 : 0)
        thread_assignments[tid] = start_idx:(start_idx + thread_size - 1)
        start_idx += thread_size
    end

    # Per-thread COO buffers
    Is = [Int[] for _ in 1:nthreads]
    Js = [Int[] for _ in 1:nthreads]
    Vs = [Float64[] for _ in 1:nthreads]

    progress  = Progress(num_emitters; dt = 1, desc = "  Bin $spectral_bin progress: ")
    completed = Threads.Atomic{Int}(0)
    inv_rays  = 1.0 / rays_per_emitter

    @threads for tid in 1:nthreads
        Il, Jl, Vl = Is[tid], Js[tid], Vs[tid]
        local_rng   = rngs[tid]
        Random.seed!(rngs[tid], seeds[tid])
        row = Dict{Int,Int}()                      # absorber → count, reused per emitter

        for global_emitter_idx in thread_assignments[tid]
            emitter_key, global_idx = all_emitters[global_emitter_idx]
            recording = rec !== nothing && global_idx in rec.ids && spectral_bin == rec.bin
            empty!(row)

            if global_idx <= num_surfaces
                coarse_index, fine_index, wall_index = emitter_key
                face = rtm.fine_mesh[coarse_index][fine_index]
                for _ in 1:rays_per_emitter
                    p_emit, dir_emit = emitSurfaceRay2D(face, wall_index, nudge, local_rng)
                    result = traceRay(rtm, p_emit, dir_emit, nudge, coarse_index, spectral_bin, local_rng)
                    result === nothing && continue
                    a = getGlobalIndex2D(surface_mapping, volume_mapping, num_surfaces, result...)
                    a == -1 && continue
                    if recording
                        push!(rec.origins[tid],   p_emit)
                        push!(rec.endpoints[tid], result[end])
                    end
                    row[a] = get(row, a, 0) + 1
                end
            else
                coarse_index, fine_index = emitter_key
                face = rtm.fine_mesh[coarse_index][fine_index]
                for _ in 1:rays_per_emitter
                    p_emit, dir_emit = emitVolumeRay2D(face, nudge, local_rng)
                    result = traceRay(rtm, p_emit, dir_emit, nudge, coarse_index, spectral_bin, local_rng)
                    result === nothing && continue
                    a = getGlobalIndex2D(surface_mapping, volume_mapping, num_surfaces, result...)
                    a == -1 && continue
                    if recording
                        push!(rec.origins[tid],   p_emit)
                        push!(rec.endpoints[tid], result[end])
                    end
                    row[a] = get(row, a, 0) + 1
                end
            end

            # Flush this emitter's row into the thread's COO buffers, already normalized to F
            for (j, c) in row
                push!(Il, global_idx); push!(Jl, j); push!(Vl, c * inv_rays)
            end

            Threads.atomic_add!(completed, 1)
            tid == 1 && update!(progress, completed[])
        end
    end
    finish!(progress)

    F_raw = sparse(reduce(vcat, Is), reduce(vcat, Js), reduce(vcat, Vs),
               num_emitters, num_emitters)
    empty!.(Is); empty!.(Js); empty!.(Vs)    # drop the live buffers
    GC.gc()                                  # now they're dead, so this actually reclaims
    return row_normalize!(F_raw, rays_per_emitter)
end

function row_normalize!(F::SparseMatrixCSC, rays_per_emitter::Int)
    rs = vec(sum(F, dims = 2))          # row sums, O(nnz)
    println("Maximum ray tracing ray loss per emitter: $(round(Int, rays_per_emitter*maximum(abs, 1 .- rs)))/$rays_per_emitter")
    rv = rowvals(F); nz = nonzeros(F)
    @inbounds for k in eachindex(nz)
        nz[k] /= rs[rv[k]]
    end
    return F
end

function group_uniform_bins(uniform_across_bin::Vector{T}; atol=1e-8, rtol=1e-8) where T
    groups = Vector{Vector{Int}}()  # each entry: indices belonging to one group
    reps   = T[]                     # representative extinction value for each group
    nonuniform = Int[]               # bins flagged -1

    for (i, v) in pairs(uniform_across_bin)
        if v < -0.1
            push!(nonuniform, i)
            continue
        end
        # find an existing group within tolerance
        idx = findfirst(r -> isapprox(r, v; atol=atol, rtol=rtol), reps)
        if idx === nothing
            push!(reps, v)
            push!(groups, [i])
        else
            push!(groups[idx], i)
        end
    end
    return groups, reps, nonuniform
end

# record ray trajectories
RayRecorder(ids::Vector{P}; bin::Integer=1, nt::P=Threads.nthreads()) where P =
    RayRecorder(ids, bin,
        [Point2{Float64}[] for _ in 1:nt],
        [Point2{Float64}[] for _ in 1:nt])

collect_rays(r::RayRecorder) =
    (reduce(vcat, r.origins), reduce(vcat, r.endpoints))