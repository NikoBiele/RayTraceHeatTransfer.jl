function exchangeRayTracing!(rtm::RayTracingDomain2D, rays_tot::P, 
                              nudge::G, verbose::Bool,
                              rec::Union{Nothing,RayRecorder}) where {G, P<:Integer}

    # Ray trace domain - returns different types based on spectral mode
    F_raw = parallelRayTracing(rtm, rays_tot, nudge, verbose; rec)

    if rtm.surfaces_only
        F_raw = F_raw[1:length(rtm.surface_mapping),1:length(rtm.surface_mapping)]
    end
    
    # Update mesh with results
    rtm.F_raw = F_raw
    rtm.F_smooth = F_raw isa Vector ?
                [Matrix{G}(undef, 2, 2) for _ in 1:length(F_raw)] : Matrix{G}(undef, 2, 2) # initialize empty
    nothing
end