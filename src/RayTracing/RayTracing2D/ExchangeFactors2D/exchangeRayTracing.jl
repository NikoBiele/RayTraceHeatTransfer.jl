function exchangeRayTracing!(rtm::RayTracingDomain2D, rays_tot::P, 
                              nudge::G, max_iters::P,
                              k_dykstra::Union{Nothing,P}, verbose::Bool,
                              rec::Union{Nothing,RayRecorder}) where {G, P<:Integer}

    # Ray trace domain - returns different types based on spectral mode
    F_raw, rays_per_emitter = parallelRayTracing(rtm, rays_tot, nudge, verbose; rec)

    if rtm.surfaces_only
        F_raw = F_raw[1:length(rtm.surface_mapping),1:length(rtm.surface_mapping)]
    end

    # Smooth exchange factors based on spectral mode
    if rtm.spectral_mode == :spectral_variable
        # Variable spectral - smooth each F matrix separately
        F_smooth = Vector{AbstractMatrix}(undef, rtm.n_spectral_bins)
        
        groups, reps, nonuniform = group_uniform_bins(rtm.uniform_across_bin)

        # loop over nonuniform
        for bin in nonuniform
            verbose && println("Smoothing F matrix for nonuniform spectral bin $bin/$(rtm.n_spectral_bins)")
            F_smooth_bin = smooth_F(
                F_raw[bin], get_w(rtm; spectral_bin=bin),
                length(rtm.surface_mapping),
                max_iters=max_iters,
                k_dykstra=k_dykstra,
                verbose=verbose,
                smooth_surfaces_only=rtm.surfaces_only
            )
            F_smooth[bin] = F_smooth_bin
        end

        for (k, idx_group) in pairs(groups)
            representative_bin = first(idx_group)
            verbose && println("Smoothing F matrix for uniform spectral bin $representative_bin/$(rtm.n_spectral_bins)")
            F_smooth_bin = smooth_F(
                F_raw[representative_bin], get_w(rtm; spectral_bin=representative_bin),
                length(rtm.surface_mapping),
                max_iters=max_iters,
                k_dykstra=k_dykstra,
                verbose=verbose,
                smooth_surfaces_only=rtm.surfaces_only
            )
            for j in idx_group
                F_smooth[j] = F_smooth_bin
            end
        end        
        
    else
        # Grey or uniform spectral - single F matrix
        if rtm.spectral_mode != :grey
            verbose && println("Smoothing single F matrix for uniform spectral extinction")
        else
            verbose && println("Smoothing single F matrix for grey extinction")
        end
        
        F_smooth = smooth_F(
            F_raw, get_w(rtm),
            length(rtm.surface_mapping),
            max_iters=max_iters,
            k_dykstra=k_dykstra,
            verbose=verbose,
            smooth_surfaces_only=rtm.surfaces_only
        )
        if rtm.spectral_mode == :spectral_variable
            F_raw = [F_raw for i in 1:rtm.n_spectral_bins]
            F_smooth = [F_smooth for i in 1:rtm.n_spectral_bins]
        end
    end

    # Update mesh with results
    rtm.F_raw = F_raw
    rtm.F_smooth = F_smooth
end