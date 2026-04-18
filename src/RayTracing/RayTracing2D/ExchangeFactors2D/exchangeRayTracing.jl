function exchangeRayTracing!(rtm::RayTracingDomain2D, rays_tot::P, tol::T, 
                              nudge::G, check_interval::P, stagnation_threshold::G, verbose::Bool) where {G, P<:Integer, T}
    
    # Ray trace domain - returns different types based on spectral mode
    F_raw, rays_per_emitter = parallelRayTracing(rtm, rays_tot, nudge, verbose)

    # Smooth exchange factors based on spectral mode
    if rtm.spectral_mode == :spectral_variable
        # Variable spectral - smooth each F matrix separately
        # println("Smoothing $rtm.n_spectral_bins F matrices for variable spectral extinction")
        # Set defaults based on the mesh's floating point precision
        F_smooth_vector = Matrix{G}[]
        
        for bin in 1:rtm.n_spectral_bins
            verbose && println("Smoothing F matrix for spectral bin $bin/$(rtm.n_spectral_bins)")
            F_smooth_bin = smoothExchangeFactors!(
                F_raw[bin], rtm, rays_per_emitter, 
                bin; max_iterations=1000, tolerance=tol,
                check_interval=check_interval,
                stagnation_threshold=stagnation_threshold,
                verbose=verbose
            )
            push!(F_smooth_vector, F_smooth_bin)
            # push!(F_smooth_uncertain_vector, F_smooth_uncertain_bin)
        end
        
        F_smooth = F_smooth_vector
        # F_smooth_uncertain = F_smooth_uncertain_vector
        
    else
        # Grey or uniform spectral - single F matrix
        if rtm.spectral_mode != :grey
            verbose && println("Smoothing single F matrix for uniform spectral extinction")
        else
            verbose && println("Smoothing single F matrix for grey extinction")
        end
        
        F_smooth = smoothExchangeFactors!(
            F_raw, rtm, rays_per_emitter, 1; max_iterations=1000,
            tolerance=tol,
            check_interval=check_interval,
            stagnation_threshold=stagnation_threshold,
            verbose=verbose
        )
    end

    # Update mesh with results
    rtm.F_raw = F_raw
    # rtm.F_raw_uncertain = F_raw_uncertain
    rtm.F_smooth = F_smooth
    # rtm.F_smooth_uncertain = F_smooth_uncertain
end