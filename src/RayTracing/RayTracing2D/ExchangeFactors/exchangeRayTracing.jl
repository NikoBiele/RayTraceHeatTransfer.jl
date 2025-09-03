function exchange_ray_tracing!(rtm::RayTracingMeshOptim, rays_tot::P, tol::G, uncertain::Bool, 
                              nudge=Float64(eps(Float32))) where {G<:AbstractFloat, P<:Integer}
    
    # ray trace domain
    F_raw, F_raw_uncertain, rays_per_emitter = parallel_ray_tracing_optimized(rtm, rays_tot, uncertain, nudge)

    # compute mappings
    surface_mapping, volume_mapping, surface_areas, volumes = get_mappings(rtm)

    # smoothen view factors with variable extinction
    F_smooth, F_smooth_uncertain = smooth_exchange_factors!(F_raw, rtm, surface_mapping, volume_mapping, rays_per_emitter, uncertain; max_iterations=100, tolerance=tol)

    # Update mesh with results
    rtm.F_raw = F_raw
    rtm.F_raw_uncertain = F_raw_uncertain
    rtm.F_smooth = F_smooth
    rtm.F_smooth_uncertain = F_smooth_uncertain
    rtm.surface_areas = surface_areas
    rtm.volumes = volumes
    rtm.surface_mapping = surface_mapping
    rtm.volume_mapping = volume_mapping
    rtm.reciprocity_satisfied = reciprocity_satisfied
    rtm.max_reciprocity_error = max_reciprocity_error
    rtm.conservation_satisfied = conservation_satisfied
    rtm.max_conservation_error = max_conservation_error
end