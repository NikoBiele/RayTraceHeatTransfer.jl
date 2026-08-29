"""
    smooth!(rtm; max_iters=1_000, k_dykstra=nothing, verbose=true, renorm=true)

Project the raw exchange factors in `rtm.F_raw` onto the set of physically
admissible exchange factors, storing the result in `rtm.F_smooth`.

Ray tracing produces exchange factors that violate reciprocity and closure to
within Monte Carlo error. `smooth!` removes those defects with `k_dykstra`
Dykstra rounds followed by alternating projections, iterating until the
reciprocity defect reaches `8 * eps(Float64)` or its noise floor.

Only `F_smooth` is written, so `F_raw` survives and the same domain may be
smoothed repeatedly under different settings.

# Keyword arguments
- `k_dykstra`: Dykstra rounds before the alternating projection finalisation.
  `0` gives pure alternating projections; larger values approach the
  minimum-norm correction. The default `nothing` picks whichever performs best
  for the problem at hand, based on the structure of `F_raw`. Rounds stop early
  once the iterate is feasible.
- `max_iters`: cap on alternating projection iterations. Reaching it warns and
  leaves `converged` false.
- `renorm`: scale the weight vector by its smallest entry before projecting.
- `verbose`: print per-bin progress and per-iteration defect bounds.

# Returns
A named tuple of diagnostics, each field a vector with one entry per spectral
bin (length 1 for grey or spectrally uniform domains):

- `converged`: whether the alternating projections reached the target
- `delta_raw`: the raw unsmoothed distance to the energy conservation/reciprocity
  feasible manifold
- `delta_smooth`: certified upper bound on the distance from `F_smooth` to the
  energy conservation/reciprocity feasible manifold
- `k_dykstra`: Dykstra rounds actually performed, which may be fewer than
  requested
- `k_ap`: alternating projection iterations performed
- `k_pcg_tot`, `k_pcg_max`: total and per-round maximum preconditioned
  conjugate gradient iterations in the dual solves; both `0` when
  `k_dykstra == 0`

For `:spectral_variable` domains, bins with identical extinction are smoothed
once and the result shared, so their entries are identical by construction.

# Example
```julia
stats = smooth!(domain)
all(stats.converged) || @warn "smoothing did not converge"
maximum(stats.delta_max) # worst-case distance to the feasible manifold across bins
```

"""
function smooth!(rtm::Union{RayTracingDomain2D,ViewFactorDomain3D};
                                 max_iters::Int=1_000,
                                 k_dykstra::Union{Nothing,Int}=nothing,
                                 verbose::Bool=true,
                                 renorm::Bool=true)

    # Smooth exchange factors based on spectral mode
    if rtm.spectral_mode == :spectral_variable && typeof(rtm) <: RayTracingDomain2D
        # Variable spectral - smooth each F matrix separately
        
        groups, reps, nonuniform = group_uniform_bins(rtm.uniform_across_bin)

        F_smooth = Vector{AbstractMatrix}(undef, rtm.n_spectral_bins)
        converged_out = Vector{Bool}(undef, rtm.n_spectral_bins)
        delta_raw_out = Vector{Float64}(undef, rtm.n_spectral_bins)
        delta_max_out = Vector{Float64}(undef, rtm.n_spectral_bins)
        k_dykstra_out = Vector{Int}(undef, rtm.n_spectral_bins)
        k_ap_out = Vector{Int}(undef, rtm.n_spectral_bins)
        k_pcg_tot_out = Vector{Int}(undef, rtm.n_spectral_bins)
        k_pcg_max_out = Vector{Int}(undef, rtm.n_spectral_bins)
        
        # loop over nonuniform
        for bin in nonuniform
            verbose && println("Smoothing F matrix for nonuniform spectral bin $bin/$(rtm.n_spectral_bins)")
            F_smooth_bin, converged_i, delta_raw_i, delta_max_i, k_dykstra_i, k_ap_i, k_pcg_tot_i, k_pcg_max_i = smooth_F(
                rtm, rtm.F_raw[bin];
                max_iters=max_iters,
                k_dykstra=k_dykstra,
                verbose=verbose,
                smooth_surfaces_only=rtm.surfaces_only,
                renorm=renorm,
                spectral_bin=bin
            )
            F_smooth[bin] = F_smooth_bin
            converged_out[bin] = converged_i
            delta_raw_out[bin] = delta_raw_i
            delta_max_out[bin] = delta_max_i
            k_dykstra_out[bin] = k_dykstra_i
            k_ap_out[bin] = k_ap_i
            k_pcg_tot_out[bin] = k_pcg_tot_i
            k_pcg_max_out[bin] = k_pcg_max_i
        end

        for (k, idx_group) in pairs(groups)
            representative_bin = first(idx_group)
            verbose && println("Smoothing F matrix for uniform spectral bin $representative_bin/$(rtm.n_spectral_bins)")
            F_smooth_bin, converged_i, delta_raw_i, delta_max_i, k_dykstra_i, k_ap_i, k_pcg_tot_i, k_pcg_max_i = smooth_F(
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
                converged_out[j] = converged_i
                delta_raw_out[j] = delta_raw_i
                delta_max_out[j] = delta_max_i
                k_dykstra_out[j] = k_dykstra_i
                k_ap_out[j] = k_ap_i
                k_pcg_tot_out[j] = k_pcg_tot_i
                k_pcg_max_out[j] = k_pcg_max_i
            end
        end        
        
    else
        # Grey or uniform spectral - single F matrix
        if rtm.spectral_mode != :grey
            verbose && println("Smoothing single F matrix for uniform spectral extinction")
        else
            verbose && println("Smoothing single F matrix for grey extinction")
        end
        
        F_smooth, converged_i, delta_raw_i, delta_max_i, k_dykstra_i, k_ap_i, k_pcg_tot_i, k_pcg_max_i = smooth_F(
            rtm, rtm.F_raw,
            max_iters=max_iters,
            k_dykstra=k_dykstra,
            verbose=verbose,
            smooth_surfaces_only=rtm.surfaces_only,
            renorm=renorm
        )
        converged_out = [converged_i]
        delta_raw_out = [delta_raw_i]
        delta_max_out = [delta_max_i]
        k_dykstra_out = [k_dykstra_i]
        k_ap_out = [k_ap_i]
        k_pcg_tot_out = [k_pcg_tot_i]
        k_pcg_max_out = [k_pcg_max_i]
    end

    rtm.F_smooth = F_smooth
    return (; converged=converged_out, delta_raw=delta_raw_out, delta_smooth=delta_max_out,
            k_dykstra=k_dykstra_out, k_ap=k_ap_out, k_pcg_tot=k_pcg_tot_out, k_pcg_max=k_pcg_max_out)
end