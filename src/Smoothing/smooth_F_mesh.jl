function smooth!(rtm::Union{RayTracingDomain2D,ViewFactorDomain3D};
                                 max_iters::Int=1_000,
                                 k_dykstra::Union{Nothing,Int}=nothing,
                                 verbose::Bool=true,
                                 renorm::Bool=true)

    # Smooth exchange factors based on spectral mode
    if rtm.spectral_mode == :spectral_variable && typeof(rtm) <: RayTracingDomain2D
        # Variable spectral - smooth each F matrix separately
        F_smooth = Vector{AbstractMatrix}(undef, rtm.n_spectral_bins)
        
        groups, reps, nonuniform = group_uniform_bins(rtm.uniform_across_bin)

        # loop over nonuniform
        for bin in nonuniform
            verbose && println("Smoothing F matrix for nonuniform spectral bin $bin/$(rtm.n_spectral_bins)")
            F_smooth_bin = smooth_F(
                rtm, rtm.F_raw[bin];
                max_iters=max_iters,
                k_dykstra=k_dykstra,
                verbose=verbose,
                smooth_surfaces_only=rtm.surfaces_only,
                renorm=renorm,
                spectral_bin=bin
            )
            F_smooth[bin] = F_smooth_bin
        end

        for (k, idx_group) in pairs(groups)
            representative_bin = first(idx_group)
            verbose && println("Smoothing F matrix for uniform spectral bin $representative_bin/$(rtm.n_spectral_bins)")
            F_smooth_bin = smooth_F(
                rtm, rtm.F_raw[representative_bin],
                max_iters=max_iters,
                k_dykstra=k_dykstra,
                verbose=verbose,
                smooth_surfaces_only=rtm.surfaces_only,
                renorm=renorm,
                spectral_bin=representative_bin
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
            rtm, rtm.F_raw,
            max_iters=max_iters,
            k_dykstra=k_dykstra,
            verbose=verbose,
            smooth_surfaces_only=rtm.surfaces_only,
            renorm=renorm
        )
    end

    rtm.F_smooth = F_smooth
    nothing
end