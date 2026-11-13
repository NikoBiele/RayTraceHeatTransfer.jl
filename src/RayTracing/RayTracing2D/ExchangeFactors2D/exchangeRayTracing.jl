function exchangeRayTracing!(rtm::RayTracingDomain2D, rays_tot::S, 
                              nudge::G, verbose::Bool,
                              rec::Union{Nothing,RayRecorder},
                              seeds::Union{UnitRange{P},Vector{P}}, rngs::Vector{<:AbstractRNG},
                              nthreads::K;
                              sampler::Symbol = :random, sobol_seed::Integer = 1,
                              ray_allocation::Symbol = :weighted
                              ) where {G, S<:Integer, P<:Integer, K<:Integer}

    # Ray trace domain - returns different types based on spectral mode
    F_raw = parallelRayTracing(rtm, rays_tot, nudge, verbose, seeds, rngs, nthreads;
                               rec, sampler, sobol_seed, ray_allocation)
                               
    # Update mesh with results
    rtm.F_raw = F_raw
    rtm.F_smooth = F_raw isa Vector ?
                [Matrix{G}(undef, 2, 2) for _ in 1:length(F_raw)] : Matrix{G}(undef, 2, 2) # initialize empty
    nothing
end