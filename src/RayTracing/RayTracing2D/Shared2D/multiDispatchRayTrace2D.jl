function (rtm::RayTracingDomain2D{VPF,VVPF,MT,VT,DIII,DII,GRID})(rays_tot::P; method::Symbol=:exchange,
                                nthreads::S=Threads.nthreads(),
                                seeds::Union{Vector{K},K,UnitRange{K},Nothing}=nothing,
                                rngs::Union{Vector{<:Random.AbstractRNG},<:Random.AbstractRNG,Nothing}=nothing,
                                nudge=nothing, verbose=true, rec=nothing) where
                                {VPF,VVPF,MT,VT,DIII,DII,P<:Integer,GRID,S<:Integer,K<:Integer}
    
    if nthreads > Threads.nthreads()
        @warn "The number of input threads is higher than the available number of the session."
    end

    if seeds === nothing
        seeds = 1:nthreads
    elseif seeds isa Integer
        seeds = seeds:(seeds + nthreads - 1)
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
        exchangeRayTracing!(rtm, rays_tot, trace_nudge, verbose, rec, seeds, rngs, nthreads)
    elseif method == :direct
        directRayTracing!(rtm, rays_tot, trace_nudge, verbose, seeds, rngs, nthreads)
    else
        error("Unknown ray tracing method: $method, must be :exchange or :direct.")
    end
end