# Record geometric ray paths once (method = :pathlength). The emitters, emission
# sampling and threading mirror computeExchangeFactorsBin; the only difference
# is that rays are walked to the wall with traceRayPath! and their cell
# sequences stored instead of one absorption tally per ray.
#
# With `chunk_rays < rays_total` the rays are recorded in chunks: each chunk is
# deposited into every bin's exchange factors and then discarded, so memory is
# bounded by the chunk size while the trace is still done once for all bins.
# The recorded paths are not kept in that case (rtm.path_store = nothing), so
# a later change of κ requires a new trace; with the default (one chunk) the
# store is kept and exchangeFactors! can be called again at no tracing cost.

function pathRayTracing!(rtm::RayTracingDomain2D, rays_total::S, nudge::G, verbose::Bool,
                         seeds::Union{UnitRange{P},Vector{P}}, rngs::Vector{<:AbstractRNG},
                         nthreads::K; chunk_rays::Integer = 10_000_000) where {S<:Integer, P<:Integer, K<:Integer, G}
    rtm.surfaces_only &&
        error("method = :pathlength requires a participating medium; use :exchange for surface-only domains")
    chunk_rays >= 1 || throw(ArgumentError("chunk_rays must be positive"))

    surface_mapping, volume_mapping, num_surfaces, num_volumes = createIndexMapping(rtm, rays_total)
    num_emitters = num_surfaces + num_volumes
    rays_per_emitter = div(rays_total, num_emitters)

    all_emitters = Vector{Tuple{Any, Int}}()
    for k in keys(surface_mapping)
        push!(all_emitters, (k, surface_mapping[k]))
    end
    for k in keys(volume_mapping)
        push!(all_emitters, (k, num_surfaces + volume_mapping[k]))
    end
    sort!(all_emitters, by = x -> x[2])

    emitters_per_thread = div(num_emitters, nthreads)
    remainder = num_emitters % nthreads
    thread_assignments = Vector{UnitRange{Int}}(undef, nthreads)
    start_idx = 1
    for tid in 1:nthreads
        thread_size = emitters_per_thread + (tid <= remainder ? 1 : 0)
        thread_assignments[tid] = start_idx:(start_idx + thread_size - 1)
        start_idx += thread_size
    end

    # chunking: rays per emitter per chunk, last chunk takes the remainder
    chunk_per_emitter = clamp(div(chunk_rays, num_emitters), 1, rays_per_emitter)
    n_chunks = cld(rays_per_emitter, chunk_per_emitter)
    chunk_sizes = fill(chunk_per_emitter, n_chunks)
    chunk_sizes[end] = rays_per_emitter - chunk_per_emitter * (n_chunks - 1)

    # per-thread buffers, reused across chunks
    cells = [Int32[] for _ in 1:nthreads]
    lens  = [Float32[] for _ in 1:nthreads]
    nseg  = [Int32[] for _ in 1:nthreads]          # segments per ray
    emit  = [Int32[] for _ in 1:nthreads]
    ends  = [Int32[] for _ in 1:nthreads]
    dirs  = [Point2{Float32}[] for _ in 1:nthreads]
    buffers = (cells, lens, nseg, emit, ends, dirs)

    # seed once; chunks continue the same streams
    for tid in 1:nthreads
        Random.seed!(rngs[tid], seeds[tid])
    end

    if n_chunks == 1
        verbose && println("Recording ray paths with $nthreads threads ($rays_per_emitter rays per emitter)")
    else
        verbose && println("Recording ray paths with $nthreads threads ($rays_per_emitter rays per emitter, ",
                           "$n_chunks chunks of ≤ $chunk_per_emitter rays per emitter)")
    end
    verbose && (progress = Progress(n_chunks * num_emitters; dt = 1, desc = "  Path recording progress: "))
    completed = Threads.Atomic{Int}(0)

    if n_chunks == 1
        store = _record_paths!(rtm, buffers, all_emitters, thread_assignments, surface_mapping, volume_mapping,
                               num_surfaces, num_volumes, rays_per_emitter, rays_per_emitter, nudge, rngs, nthreads,
                               verbose, verbose ? progress : nothing, completed)
        verbose && finish!(progress)
        empty!.(cells); empty!.(lens); empty!.(nseg); empty!.(emit); empty!.(ends); empty!.(dirs)
        GC.gc()

        lost = num_emitters * rays_per_emitter - n_rays(store)
        verbose && println("Recorded $(n_rays(store)) rays, $(n_segments(store)) segments ",
                           "($(round(8 * n_segments(store) / 2^20; digits = 1)) MiB); lost rays: $lost")

        rtm.path_store = store
        exchangeFactors!(rtm; verbose)
        return nothing
    end

    # ---- chunked: deposit each chunk into every bin, keep only the triplets ----
    N = num_surfaces + num_volumes
    n_bins = rtm.n_spectral_bins
    βT = permutedims(_extinction_table(rtm, num_volumes, n_bins))   # bin-major for the inner loop
    acc = [spzeros(Float64, Int, N, N) for _ in 1:n_bins]   # unnormalised row sums, merged per chunk
    recorded = 0; segments = 0

    for (ic, rays_this_chunk) in enumerate(chunk_sizes)
        store = _record_paths!(rtm, buffers, all_emitters, thread_assignments, surface_mapping, volume_mapping,
                               num_surfaces, num_volumes, rays_this_chunk, rays_per_emitter, nudge, rngs, nthreads,
                               verbose, verbose ? progress : nothing, completed)
        empty!.(cells); empty!.(lens); empty!.(nseg); empty!.(emit); empty!.(ends); empty!.(dirs)

        recorded += n_rays(store); segments += n_segments(store)

        I, J, V = _deposit_all_bins(store, βT, num_surfaces, N)
        Threads.@threads for k in 1:n_bins
            acc[k] = acc[k] + sparse(I[k], J[k], V[k], N, N)   # bounded by nnz(F), not by chunk count
        end
        store = nothing
        GC.gc()
    end
    verbose && finish!(progress)

    lost = num_emitters * rays_per_emitter - recorded
    verbose && println("Recorded $recorded rays in $n_chunks chunks, $segments segments ",
                       "(≤ $(round(8 * div(segments, n_chunks) / 2^20; digits = 1)) MiB held at once); lost rays: $lost")

    F_bins = Vector{SparseMatrixCSC{Float64,Int}}(undef, n_bins)
    Threads.@threads for k in 1:n_bins
        I, J, V = findnz(acc[k])
        F_bins[k] = _assemble_bin(I, J, V, N)
    end
    verbose && println("Exchange factors from $recorded recorded rays for $n_bins bin(s)")

    rtm.path_store = nothing                      # nothing kept: re-binning needs a new trace
    _set_F_raw!(rtm, F_bins)
    return nothing
end

# Record `rays_this_chunk` rays from every emitter into the per-thread buffers
# and return them as a RayPathStore. Buffers must be empty on entry; the rngs
# are used as they are (no reseeding), so successive calls continue the streams.
function _record_paths!(rtm::RayTracingDomain2D, buffers, all_emitters, thread_assignments,
                        surface_mapping, volume_mapping, num_surfaces::Int, num_volumes::Int,
                        rays_this_chunk::Int, rays_per_emitter::Int, nudge, rngs, nthreads::Int,
                        verbose::Bool, progress, completed)
    cells, lens, nseg, emit, ends, dirs = buffers

    @threads for tid in 1:nthreads
        cl, ll, nl, el, wl, dl = cells[tid], lens[tid], nseg[tid], emit[tid], ends[tid], dirs[tid]
        local_rng = rngs[tid]

        for global_emitter_idx in thread_assignments[tid]
            emitter_key, global_idx = all_emitters[global_emitter_idx]
            is_surface = global_idx <= num_surfaces
            coarse_index::Int = emitter_key[1]
            fine_index::Int   = emitter_key[2]
            wall_index::Int   = is_surface ? emitter_key[3] : 0
            face = rtm.fine_mesh[coarse_index][fine_index]

            for _ in 1:rays_this_chunk
                p_emit, dir_emit = is_surface ? emitSurfaceRay2D(face, wall_index, nudge, local_rng) :
                                                emitVolumeRay2D(face, nudge, local_rng)
                n_before = length(cl)
                result = traceRayPath!(cl, ll, rtm, p_emit, dir_emit, nudge, coarse_index,
                                       volume_mapping, num_surfaces)
                a = result === 0 ? -1 :
                    getGlobalIndex2D(surface_mapping, volume_mapping, num_surfaces, result[1], result[2], result[3], result[4])
                if a == -1
                    resize!(cl, n_before); resize!(ll, n_before)   # lost ray: discard its segments
                    continue
                end
                push!(nl, Int32(length(cl) - n_before))
                push!(el, Int32(global_idx))
                push!(wl, Int32(a))
                push!(dl, Point2{Float32}(dir_emit[1], dir_emit[2]))
            end

            verbose && Threads.atomic_add!(completed, 1)
            verbose && (tid == 1 && update!(progress, completed[]))
        end
    end

    seg_cell    = reduce(vcat, cells)
    seg_len     = reduce(vcat, lens)
    ray_nseg    = reduce(vcat, nseg)
    ray_emitter = reduce(vcat, emit)
    ray_end     = reduce(vcat, ends)
    ray_dir     = reduce(vcat, dirs)
    ray_start   = Vector{Int}(undef, length(ray_nseg) + 1)
    ray_start[1] = 1
    cumsum!(view(ray_start, 2:length(ray_start)), Int.(ray_nseg))
    ray_start[2:end] .+= 1

    return RayPathStore(seg_cell, seg_len, ray_start, ray_emitter, ray_end, ray_dir,
                        num_surfaces, num_volumes, rays_per_emitter)
end

n_rays(s::RayPathStore) = length(s.ray_emitter)
n_segments(s::RayPathStore) = length(s.seg_cell)

# β per global volume index, one column per bin
function _extinction_table(rtm::RayTracingDomain2D, nv::Int, n_bins::Int)
    β = zeros(nv, n_bins)
    for ((ci, fi), v) in rtm.volume_mapping
        face = rtm.fine_mesh[ci][fi]
        for k in 1:n_bins
            β[v, k] = face.kappa_g isa Vector ? face.kappa_g[k] + face.sigma_s_g[k] :
                                                face.kappa_g + face.sigma_s_g
        end
    end
    return β
end

function _set_F_raw!(rtm::RayTracingDomain2D, F_bins::Vector{SparseMatrixCSC{Float64,Int}})
    is_spectral = rtm.spectral_mode != :grey
    F_raw = is_spectral ? Vector{AbstractMatrix}(F_bins) : F_bins[1]
    rtm.F_raw = F_raw
    rtm.F_smooth = F_raw isa Vector ?
        [Matrix{Float64}(undef, 2, 2) for _ in 1:length(F_raw)] : Matrix{Float64}(undef, 2, 2)
    return F_raw
end

"""
    exchangeFactors!(rtm::RayTracingDomain2D; verbose = true)

Compute `rtm.F_raw` for every spectral bin from the recorded ray paths
(`rtm.path_store`) and the current extinction coefficients of the fine
faces. Along each ray the fraction absorbed in cell c is
e^{-τ}·(1 − e^{-β_c ℓ_c}) with τ the optical depth accumulated before it,
and what remains reaches the wall the ray ended on; every ray deposits
exactly 1, so rows sum to 1 up to rounding. Changing κ and calling this
again costs one pass over the store, not a new trace.
"""
function exchangeFactors!(rtm::RayTracingDomain2D; verbose::Bool = true)
    store = rtm.path_store
    store === nothing && error("no recorded ray paths: call domain(N_rays; method = :pathlength) first " *
                               "(a chunked trace, chunk_rays < N_rays, does not keep the paths)")

    ns, nv = store.num_surfaces, store.num_volumes
    N = ns + nv
    n_bins = rtm.n_spectral_bins
    βT = permutedims(_extinction_table(rtm, nv, n_bins))

    I, J, V = _deposit_all_bins(store, βT, ns, N)
    F_bins = Vector{SparseMatrixCSC{Float64,Int}}(undef, n_bins)
    Threads.@threads for k in 1:n_bins
        F_bins[k] = _assemble_bin(I[k], J[k], V[k], N)
    end
    verbose && println("Exchange factors from $(n_rays(store)) recorded rays for $n_bins bin(s)")

    return _set_F_raw!(rtm, F_bins)
end

# Normalise each row by its emitter's kept ray count and assemble; sparse()
# sums duplicate triplets, which is what merges the chunks.
function _assemble_bin(I::Vector{Int}, J::Vector{Int}, V::Vector{Float64}, N::Int)
    rowsum = zeros(N)
    @inbounds for i in eachindex(V)
        rowsum[I[i]] += V[i]
    end
    @inbounds for i in eachindex(V)
        V[i] /= rowsum[I[i]]
    end
    return sparse(I, J, V, N, N)
end

# Deposit one store into every bin in a single pass over the segments. Rays
# are stored in emitter order, so emitters are split into contiguous blocks,
# one per thread, each with its own dense row buffer (n_bins × N) that is
# flushed to per-bin triplet lists after each emitter. Along a ray the
# transmission t_k = e^{-τ_k} is carried per bin and updated multiplicatively,
# so each segment costs one exp per bin. Returns per-bin (I, J, V) with
# UNNORMALISED row sums; duplicate (I, J) pairs across blocks/chunks are summed
# by sparse() in _assemble_bin. βT is n_bins × n_volumes (bin-major).
function _deposit_all_bins(store::RayPathStore, βT::Matrix{Float64}, ns::Int, N::Int)
    n_bins = size(βT, 1)
    seg_cell, seg_len, ray_start = store.seg_cell, store.seg_len, store.ray_start
    ray_emitter, ray_end = store.ray_emitter, store.ray_end

    # emitter e owns rays e_start[e]:e_start[e+1]-1
    cnt = zeros(Int, N)
    for e in ray_emitter
        cnt[e] += 1
    end
    e_start = cumsum(vcat(1, cnt))

    nt = Threads.nthreads()
    blocks = collect(Iterators.partition(1:N, cld(N, nt)))
    nb = length(blocks)
    Ib = [[Int[] for _ in 1:n_bins] for _ in 1:nb]
    Jb = [[Int[] for _ in 1:n_bins] for _ in 1:nb]
    Vb = [[Float64[] for _ in 1:n_bins] for _ in 1:nb]

    Threads.@threads for b in 1:nb
        row     = zeros(n_bins, N)
        t       = ones(n_bins)
        marker  = falses(N)
        touched = Int[]
        Is, Js, Vs = Ib[b], Jb[b], Vb[b]

        @inbounds for e in blocks[b]
            for r in e_start[e]:e_start[e+1]-1
                fill!(t, 1.0)
                for s in ray_start[r]:ray_start[r+1]-1
                    c = Int(seg_cell[s])
                    v = c - ns
                    ℓ = Float64(seg_len[s])
                    if !marker[c]
                        marker[c] = true; push!(touched, c)
                    end
                    @simd for k in 1:n_bins
                        ek = exp(-βT[k, v] * ℓ)
                        row[k, c] += t[k] * (1.0 - ek)
                        t[k] *= ek
                    end
                end
                w = Int(ray_end[r])
                if !marker[w]
                    marker[w] = true; push!(touched, w)
                end
                @simd for k in 1:n_bins
                    row[k, w] += t[k]
                end
            end
            for c in touched
                for k in 1:n_bins
                    val = row[k, c]
                    if val != 0.0
                        push!(Is[k], e); push!(Js[k], c); push!(Vs[k], val)
                        row[k, c] = 0.0
                    end
                end
                marker[c] = false
            end
            empty!(touched)
        end
    end

    I = [reduce(vcat, (Ib[b][k] for b in 1:nb)) for k in 1:n_bins]
    J = [reduce(vcat, (Jb[b][k] for b in 1:nb)) for k in 1:n_bins]
    V = [reduce(vcat, (Vb[b][k] for b in 1:nb)) for k in 1:n_bins]
    return I, J, V
end