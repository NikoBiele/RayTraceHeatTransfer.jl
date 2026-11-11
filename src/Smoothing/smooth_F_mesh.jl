"""
    smooth!(rtm; k_ap=nothing, k_dykstra=1, verbose=true, renorm=true)

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
  minimum-norm correction. The value is only used where the structure of
  `F_raw` leaves a choice. Input traced with `method = :pathlength` is
  smoothed by alternating projections alone. So is optically thick 2D input,
  whose gas escape fraction (the fraction of gas emission absorbed at the
  walls) lies below a threshold that falls slowly with the number of elements,
  and sparse input too large to densify; a dense input whose a priori rates
  make Dykstra cheaper per decade of defect reduction than alternating
  projections (near-parallel-plate enclosures) runs Dykstra to its floor
  regardless of `k_dykstra`. In between, the prescribed rounds are run;
  the default of 1 gives the orthogonal-projection start, and more rounds
  move the result toward the minimum-norm correction. Rounds stop early once
  the iterate is feasible or its defect stops improving.
- `renorm`: scale the weight vector by its smallest entry before projecting.
- `verbose`: print per-bin progress and per-iteration defect bounds.
- `k_ap`: alternating projection iterations. With the default `nothing` they
  run until the iterate converges, reaches its floor, or stops improving (the
  best defect falls by less than 0.1 % over 1 000 iterations), with a ceiling
  of 10^6. An integer caps them instead. Stopping at a cap or a stall warns and
  leaves `converged` false. The same applies to a domain with a
  `directional_model`. When the a priori rate predicts more than 10 000
  iterations, `smooth!` warns before it starts.

With a `directional_model` set, the angular exchange factors `G_raw` are smoothed
instead (per-bin emission shares, reversed-bin reciprocity, zero pattern kept),
independently for every pair of opposite direction bins; `G_smooth` is written and
`F_smooth` is their sum over direction bins. Only alternating projections are used
there, so `k_dykstra` is ignored and reported as 0.

# Returns
A named tuple of diagnostics, each field a vector with one entry per spectral
bin (length 1 for grey or spectrally uniform domains):

- `converged`: whether the alternating projections reached the target
- `delta_raw`: the raw unsmoothed distance to the energy conservation/reciprocity
  feasible manifold
- `delta_smooth`: distance from `F_smooth` to the energy conservation/reciprocity
  feasible manifold: a certified upper bound when alternating projections
  finalised, the exact dual distance when Dykstra reached its floor on its own
- `k_dykstra`: Dykstra rounds actually performed, which may differ from the
  value requested (see `k_dykstra` above)
- `k_ap`: alternating projection iterations performed; `0` when Dykstra reached
  its floor and no finalisation was needed
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
function smooth!(rtm::Union{RayTracingDomain2D,SurfaceDomain3D};
                                 k_ap::Union{Nothing,Int}=nothing,
                                 k_dykstra::Int=1,
                                 verbose::Bool=true,
                                 renorm::Bool=true,
                                 keep_F_raw::Bool=true)

    # directional domain: smooth the angular exchange factors, F_smooth = Σₐ G_smooth[a];
    # k_ap = nothing (the default) stops alternating projections on progress, an Int caps them
    if rtm isa RayTracingDomain2D && rtm.directional_model !== nothing
        return smooth_G!(rtm; k_ap = k_ap, verbose = verbose, keep_F_raw = keep_F_raw)
    end

    # A priori alternating-projection rate, one value per spectral bin, from the
    # blind macro-surfaces (straight walls in 2D, planar faces in 3D). In 2D it
    # also feeds the cost rule that sets the Dykstra budget; in 3D it only informs
    # the output below, and every bin gets `nothing` for the routing.
    rho_bins = nothing
    pred = predict_rho_data(rtm)
    if rtm isa RayTracingDomain2D
        rho_bins = pred.rho_bins
    end
    if verbose
        if pred.exact
            println("A priori alternating-projection rate ρ_AP ≈ $(round(pred.rho_est; sigdigits=6)) (worst bin)")
        else
            println("A priori alternating-projection rate ρ_AP ≥ $(round(pred.rho_est; sigdigits=6)) (worst bin, blind-set search truncated)")
        end
    end

    # Warn when alternating projections are predicted to need many iterations:
    # from the raw defect to the floor at the predicted rate. A first-order
    # estimate, exact near the parallel-plate limit, that fine meshes away from
    # it can exceed several-fold.
    if rtm.F_raw isa Vector
        N_el = size(rtm.F_raw[1], 1)
    else
        N_el = size(rtm.F_raw, 1)
    end
    guard = sqrt(N_el) * 8 * eps(Float64)
    k_slow = Int[]
    n_unbounded = 0
    for rho in pred.rho_bins
        if rho >= 1
            n_unbounded += 1
        else
            k_est = ceil(Int, log(guard / max(pred.delta_init, guard)) / log(rho))
            if k_est > 10_000
                push!(k_slow, k_est)
            end
        end
    end
    if !isempty(k_slow) || n_unbounded > 0
        msg = "Alternating projections are predicted to converge slowly (ρ_AP up to $(round(pred.rho_est; sigdigits=6))):"
        if !isempty(k_slow)
            msg *= "\n  about $(maximum(k_slow)) iterations in the slowest of $(length(k_slow)) bin(s) above 10 000;" *
                   " an estimate that fine meshes can exceed several-fold."
        end
        if n_unbounded > 0
            msg *= "\n  $(n_unbounded) bin(s) without a finite estimate (exactly bipartite blind set)."
        end
        msg *= "\n  With k_ap = nothing (the default) smoothing runs until it converges or stalls; an explicit k_ap caps it."
        @warn msg
    end
    rho_for(bin) = rho_bins === nothing ? nothing : rho_bins[bin]

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
                k_ap=k_ap,
                k_dykstra=k_dykstra,
                verbose=verbose,
                smooth_surfaces_only=rtm.surfaces_only,
                renorm=renorm,
                spectral_bin=bin,
                rho_ap=rho_for(bin)
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
                k_ap=k_ap,
                k_dykstra=k_dykstra,
                verbose=verbose,
                smooth_surfaces_only=rtm.surfaces_only,
                renorm=renorm,
                spectral_bin=representative_bin,
                rho_ap=rho_for(representative_bin)
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
            k_ap=k_ap,
            k_dykstra=k_dykstra,
            verbose=verbose,
            smooth_surfaces_only=rtm.surfaces_only,
            renorm=renorm,
            rho_ap=rho_for(1)
        )
        converged_out = [converged_i]
        delta_raw_out = [delta_raw_i]
        delta_max_out = [delta_max_i]
        k_dykstra_out = [k_dykstra_i]
        k_ap_out = [k_ap_i]
        k_pcg_tot_out = [k_pcg_tot_i]
        k_pcg_max_out = [k_pcg_max_i]
    end

    keep_F_raw || (rtm.F_raw = Matrix{Float64}(undef, 2, 2))
    rtm.F_smooth = F_smooth
    return (; converged=converged_out, delta_raw=delta_raw_out, delta_smooth=delta_max_out,
            k_dykstra=k_dykstra_out, k_ap=k_ap_out, k_pcg_tot=k_pcg_tot_out, k_pcg_max=k_pcg_max_out)
end