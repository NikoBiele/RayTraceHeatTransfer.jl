function exchangeRayTracing!(rtm::RayTracingDomain2D, rays_tot::P, 
                              nudge::G, verbose::Bool,
                              rec::Union{Nothing,RayRecorder}) where {G, P<:Integer}

    # Ray trace domain - returns different types based on spectral mode
    F_raw = parallelRayTracing(rtm, rays_tot, nudge, verbose; rec)
    
    # Update mesh with results
    rtm.F_raw = F_raw
    rtm.F_smooth = F_raw isa Vector ?
                [Matrix{G}(undef, 2, 2) for _ in 1:length(F_raw)] : Matrix{G}(undef, 2, 2) # initialize empty
    pred = predict_rho_data(rtm)
    if pred.rho_est > 0.99 && pred.exact
        N     = size(F_raw isa Vector ? F_raw[1] : F_raw, 1)
        guard = sqrt(N) * 8 * eps(Float64)
        iters = trunc(Int, log(guard / pred.delta_init) / log(pred.rho_est))
        @warn "Parallel-plate-like coupling detected in the traced exchange factors.\n" *
              "Pure alternating-projection smoothing predicted to converge at ρ ≈ $(round(pred.rho_est; digits=8)).\n" *
              "Consider pure Dykstra smoothing if feasible: 'smooth!(mesh; k_dykstra=1000, max_iters=0)',\n" *
              "or raise 'max_iters' to at least $(5*iters) for AP to reach the target."
    end
    if pred.rho_est > 0.99 && !pred.exact
        @warn "Parallel-plate-like coupling detected in the traced exchange factors.\n"
                "ρ ≥ $(round(pred.rho_est; digits=6)); slow alternating-projection convergence possible.\n" *
                "The estimate is a lower bound; if alternating-projection smoothing stalls, rerun with more rays,\n"
                "or use Dykstra, if feasible: 'smooth!(mesh; k_dykstra=1000, max_iters=0)'."
    end
    nothing
end