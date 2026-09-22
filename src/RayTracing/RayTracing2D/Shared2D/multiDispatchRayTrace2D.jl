function (rtm::RayTracingDomain2D{VPF,VVPF,MT,VT,DIII,DII,GRID})(rays_tot::P; method::Symbol=:exchange,
                                nthreads::S=Threads.nthreads(),
                                seeds::Union{Vector{K},K,UnitRange{K},Nothing}=nothing,
                                rngs::Union{Vector{<:Random.AbstractRNG},<:Random.AbstractRNG,Nothing}=nothing,
                                nudge=nothing, verbose=true, rec=nothing, chunk_rays::Integer=10_000_000,
                                sampler::Union{Symbol,Nothing}=nothing) where
                                {VPF,VVPF,MT,VT,DIII,DII,P<:Integer,GRID,S<:Integer,K<:Integer}
    
    if nthreads > Threads.nthreads()
        @warn "The number of input threads is higher than the available number of the session."
    end

    # Sampler: Sobol points wherever a ray consumes a fixed number of random numbers,
    # pseudorandom numbers otherwise. `sampler = nothing` picks this automatically.
    sobol_methods = (:exchange, :pathlength,)
    if sampler === nothing
        sampler = method in sobol_methods ? :sobol : :random
    end
    sampler in (:sobol, :random) ||
        error("Unknown sampler: $sampler, must be :sobol or :random.")
    sobol_seed = 1
    if sampler == :sobol
        method in sobol_methods ||
            error("sampler = :sobol is not available for method = :$method; it requires a tracer in which " *
                  "every ray consumes a fixed number of random numbers ($(join(sobol_methods, ", "))). " *
                  "Omit `sampler` or pass sampler = :random.")
        rngs === nothing ||
            error("`rngs` has no meaning with Sobol sampling (the default for method = :$method): the points " *
                  "come from one Sobol sequence per emitter. Pass sampler = :random to use your own generators.")
        seeds === nothing || seeds isa Integer ||
            error("Sobol sampling (the default for method = :$method) takes a single integer seed, which selects " *
                  "the realisation; the result does not depend on the number of threads. " *
                  "Pass sampler = :random to use one seed per thread.")
        if seeds !== nothing
            seeds >= 1 || error("an integer seed must be ≥ 1, got $seeds")
            sobol_seed = Int(seeds)
        end
        seeds = nothing      # the per-thread seeds below are unused by the Sobol sampler
    end

    if seeds === nothing
        seeds = 1:nthreads
    elseif seeds isa Integer
        # the seeds-th block of nthreads seeds, so different integers never share a thread stream
        seeds >= 1 || error("an integer seed must be ≥ 1, got $seeds")
        seeds = ((seeds - 1) * nthreads + 1):(seeds * nthreads)
    end
    length(seeds) == nthreads ||
        error("got $(length(seeds)) seeds for $nthreads threads; supply one per thread, a single starting seed, or none")
    length(unique(seeds)) == nthreads ||
        error("seeds must be distinct; threads sharing a seed trace identical rays")

    if rngs === nothing
        rngs = [Random.Xoshiro(s) for s in seeds]
    elseif rngs isa Random.AbstractRNG
        rngs = [deepcopy(rngs) for _ in 1:nthreads]
    end
    length(rngs) == nthreads ||
        error("got $(length(rngs)) rngs for $nthreads threads")
    length(unique(objectid.(rngs))) == nthreads ||
        error("rngs must be distinct objects; `fill(Xoshiro(1), n)` aliases one RNG across all threads")
        
    # Extract floating point type from the mesh vertices (Point2{G} where G is the float type)
    FloatType = eltype(rtm.fine_mesh[1][1].T_in_g) # Gets G from PolyVolume2D{G}
    
    # Set defaults based on the mesh's floating point precision
    trace_nudge = nudge === nothing ? 10_000 * eps(FloatType) : nudge

    if method == :exchange
        exchangeRayTracing!(rtm, rays_tot, trace_nudge, verbose, rec, seeds, rngs, nthreads;
                            sampler=sampler, sobol_seed=sobol_seed)
    elseif method == :direct
        directRayTracing!(rtm, rays_tot, trace_nudge, verbose, seeds, rngs, nthreads)
    elseif method == :pathlength
        pathRayTracing!(rtm, rays_tot, trace_nudge, verbose, seeds, rngs, nthreads;
                        chunk_rays=chunk_rays, sampler=sampler, sobol_seed=sobol_seed)
    else
        error("Unknown ray tracing method: $method, must be :exchange, :direct or :pathlength.")
    end
end