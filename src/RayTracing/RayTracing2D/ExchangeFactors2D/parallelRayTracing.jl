function parallelRayTracing(rtm::RayTracingDomain2D, rays_total::S, 
                                        nudge::G, verbose::Bool,
                                        seeds::Union{UnitRange{P},Vector{P}}, rngs::Vector{<:AbstractRNG},
                                        nthreads::K;
                                        rec=nothing, sampler::Symbol=:random,
                                        sobol_seed::Integer=1,
                                        ray_allocation::Symbol=:weighted) where {S<:Integer, P<:Integer, K<:Integer, G}

    surface_mapping, volume_mapping, num_surfaces, num_volumes = createIndexMapping(rtm, rays_total)
    num_emitters = num_surfaces + num_volumes

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
            ray_counts = emitter_ray_counts(rtm, rays_total, num_emitters, ray_allocation, bin)
            F_raw_bin = computeExchangeFactorsBin(
                rtm, ray_counts, nudge, bin,
                surface_mapping, volume_mapping, num_surfaces,
                num_volumes, num_emitters, verbose, rec, seeds, rngs, nthreads;
                sampler, sobol_seed
            )
            F_raw_vector[bin] = F_raw_bin
        end

        # loop over groups
        for (k, idx_group) in pairs(groups)
            representative_bin = first(idx_group)
            verbose && println("Computing F matrix for uniform spectral bin $representative_bin/$n_bins")
            ray_counts = emitter_ray_counts(rtm, rays_total, num_emitters, ray_allocation, representative_bin)
            F_raw_bin = computeExchangeFactorsBin(
                rtm, ray_counts, nudge, representative_bin,
                surface_mapping, volume_mapping, num_surfaces,
                num_volumes, num_emitters, verbose, rec, seeds, rngs, nthreads;
                sampler, sobol_seed
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
        
        ray_counts = emitter_ray_counts(rtm, rays_total, num_emitters, ray_allocation, 1)
        F_raw = computeExchangeFactorsBin(
            rtm, ray_counts, nudge, 1,  # Use bin 1 (doesn't matter for uniform)
            surface_mapping, volume_mapping, num_surfaces,
            num_volumes, num_emitters, verbose, rec, seeds, rngs, nthreads;
            sampler, sobol_seed
        )
        
        if rtm.surfaces_only
            return F_raw[1:length(rtm.surface_mapping),1:length(rtm.surface_mapping)]
        else
            return F_raw
        end
    end
end

function computeExchangeFactorsBin(rtm::RayTracingDomain2D, ray_counts::Vector{Int},
                                 nudge::G, spectral_bin::P,
                                 surface_mapping, volume_mapping, num_surfaces,
                                 num_volumes, num_emitters, verbose, rec,
                                 seeds::Union{UnitRange{K},Vector{K}}, rngs::Vector{<:AbstractRNG},
                                 nthreads::P; sampler::Symbol = :random,
                                 sobol_seed::Integer = 1) where {P<:Integer, K<:Integer, G}

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

    # Partition emitters across threads — disjoint ranges with about equal numbers of
    # rays, so each thread is independent
    thread_assignments = partition_emitters(ray_counts, nthreads)

    # Per-thread COO buffers
    Is = [Int[] for _ in 1:nthreads]
    Js = [Int[] for _ in 1:nthreads]
    Vs = [Float64[] for _ in 1:nthreads]

    verbose && (progress  = Progress(num_emitters; dt = 1, desc = "  Bin $spectral_bin progress: "))
    verbose && (completed = Threads.Atomic{Int}(0))

    # random numbers: one Sobol source per emitter (:sobol) or one stream per thread (:random).
    per_emitter = sampler == :sobol
    ray_rngs = per_emitter ? sobolEmitterRNGs2D(rtm, all_emitters, num_surfaces, true, sobol_seed, spectral_bin) : rngs

    @threads for tid in 1:nthreads
        Il, Jl, Vl = Is[tid], Js[tid], Vs[tid]
        per_emitter || Random.seed!(rngs[tid], seeds[tid])
        row = Dict{Int,Int}()                      # absorber → count, reused per emitter

        for global_emitter_idx in thread_assignments[tid]
            emitter_key, global_idx = all_emitters[global_emitter_idx]
            local_rng = per_emitter ? ray_rngs[global_emitter_idx] : ray_rngs[tid]   # the emitter's Sobol source, or the thread's stream
            recording = rec !== nothing && global_idx in rec.ids && spectral_bin == rec.bin
            empty!(row)

            is_surface = global_idx <= num_surfaces
            coarse_index::Int = emitter_key[1]
            fine_index::Int   = emitter_key[2]
            wall_index::Int   = is_surface ? emitter_key[3] : 0
            face = rtm.fine_mesh[coarse_index][fine_index]

            for _ in 1:ray_counts[global_idx]
                nextRay!(local_rng)                        # next Sobol point (no-op for ordinary generators)
                p_emit, dir_emit = is_surface ? emitSurfaceRay2D(face, wall_index, nudge, local_rng) :
                                                emitVolumeRay2D(face, nudge, local_rng)
                result = traceRay(rtm, p_emit, dir_emit, nudge, coarse_index, spectral_bin, local_rng)
                result === 0 && continue
                a = getGlobalIndex2D(surface_mapping, volume_mapping, num_surfaces, result[1], result[2], result[3], result[4])
                a == -1 && continue
                if recording
                    push!(rec.origins[tid],   p_emit)
                    push!(rec.endpoints[tid], result[end])
                end
                row[a] = get(row, a, 0) + 1
            end

            # Flush this emitter's row into the thread's COO buffers, already normalized to F
            inv_rays = 1.0 / ray_counts[global_idx]
            for (j, c) in row
                push!(Il, global_idx); push!(Jl, j); push!(Vl, c * inv_rays)
            end

            verbose && Threads.atomic_add!(completed, 1)
            verbose && (tid == 1 && update!(progress, completed[]))
        end
    end
    verbose && finish!(progress)

    F_raw = sparse(reduce(vcat, Is), reduce(vcat, Js), reduce(vcat, Vs),
               num_emitters, num_emitters)
    empty!.(Is); empty!.(Js); empty!.(Vs)    # drop the live buffers
    GC.gc()                                  # now they're dead, so this actually reclaims
    return row_normalize!(F_raw, ray_counts, verbose)
end

# Rays per emitter for one spectral bin. :uniform gives every emitter the same
# number. :weighted gives rays in proportion to the reciprocity weight w of the
# bin (area for walls, 4 beta V for gas cells), with the gas share of the rays
# capped at its equal-count share N_g/N: the weighting may move rays from the gas
# to the walls, which helps in optically thin media, but never from the walls
# into the gas, which would starve the wall rows in optically thick media. Within
# the walls and within the gas the rays stay in proportion to w. Every emitter
# gets at least one ray.
function emitter_ray_counts(rtm::RayTracingDomain2D, rays_total::Integer, num_emitters::Int,
                            ray_allocation::Symbol, spectral_bin::Integer)
    if ray_allocation == :uniform
        return fill(div(rays_total, num_emitters), num_emitters)
    end
    w = get_w(rtm; spectral_bin = spectral_bin)
    length(w) == num_emitters ||
        error("got $(length(w)) weights for $num_emitters emitters")
    n_surf = length(rtm.surface_mapping)
    gas = (n_surf + 1):num_emitters
    if rtm.surfaces_only
        # the gas rows are dropped after the trace, so the gas cells get one ray each
        w[gas] .= 0.0
    end
    w_total = sum(w)
    w_total > 0 || error("all emitter weights are zero in spectral bin $spectral_bin")
    w_gas = sum(w[gas])
    share_count = length(gas) / num_emitters
    if w_gas / w_total > share_count
        # cap the gas at its equal-count share, keeping w within each class
        w[1:n_surf] .*= (1.0 - share_count) / (w_total - w_gas)
        w[gas] .*= share_count / w_gas
        w_total = sum(w)
    end
    counts = Vector{Int}(undef, num_emitters)
    for i in 1:num_emitters
        counts[i] = max(1, round(Int, rays_total * w[i] / w_total))
    end
    return counts
end

# Contiguous emitter ranges, one per thread, with about equal numbers of rays. An
# emitter goes to the current thread while the midpoint of its rays lies within
# that thread's share. Equal counts reproduce the plain split by emitter number.
function partition_emitters(ray_counts::Vector{Int}, nthreads::Integer)
    num_emitters = length(ray_counts)
    ranges = Vector{UnitRange{Int}}(undef, nthreads)
    if all(==(ray_counts[1]), ray_counts)
        emitters_per_thread = div(num_emitters, nthreads)
        remainder = num_emitters % nthreads
        start_idx = 1
        for tid in 1:nthreads
            thread_size = emitters_per_thread
            if tid <= remainder
                thread_size += 1
            end
            ranges[tid] = start_idx:(start_idx + thread_size - 1)
            start_idx += thread_size
        end
        return ranges
    end
    total = sum(ray_counts)
    start_idx = 1
    cum = 0
    for tid in 1:nthreads
        stop_idx = start_idx - 1
        if tid == nthreads
            stop_idx = num_emitters
        else
            share = total * tid / nthreads
            while stop_idx < num_emitters && cum + ray_counts[stop_idx + 1] / 2 <= share
                stop_idx += 1
                cum += ray_counts[stop_idx]
            end
        end
        ranges[tid] = start_idx:stop_idx
        start_idx = stop_idx + 1
    end
    return ranges
end

function row_normalize!(F::SparseMatrixCSC, rays_per_emitter::Int, verbose::Bool)
    rs = vec(sum(F, dims = 2))          # row sums, O(nnz)
    verbose && println("Maximum ray tracing ray loss per emitter: $(round(Int, rays_per_emitter*maximum(abs, 1 .- rs)))/$rays_per_emitter")
    rv = rowvals(F); nz = nonzeros(F)
    @inbounds for k in eachindex(nz)
        nz[k] /= rs[rv[k]]
    end
    return F
end

# Row normalisation for emitters with individual ray counts; the printed loss is
# that of the emitter losing the largest fraction of its rays
function row_normalize!(F::SparseMatrixCSC, ray_counts::Vector{Int}, verbose::Bool)
    rs = vec(sum(F, dims = 2))          # row sums, O(nnz)
    if verbose
        worst = argmax(abs.(1 .- rs))
        lost = round(Int, ray_counts[worst] * abs(1 - rs[worst]))
        println("Maximum ray tracing ray loss per emitter: $lost/$(ray_counts[worst])")
    end
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