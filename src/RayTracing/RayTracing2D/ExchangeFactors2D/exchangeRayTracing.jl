function exchangeRayTracing!(rtm::RayTracingDomain2D, rays_tot::S, 
                              nudge::G, verbose::Bool,
                              rec::Union{Nothing,RayRecorder},
                              seeds::Union{UnitRange{P},Vector{P}}, rngs::Vector{<:AbstractRNG},
                              nthreads::K;
                              sampler::Symbol = :random, sobol_seed::Integer = 1
                              ) where {G, S<:Integer, P<:Integer, K<:Integer}

    # Ray trace domain - returns different types based on spectral mode
    F_raw = parallelRayTracing(rtm, rays_tot, nudge, verbose, seeds, rngs, nthreads;
                               rec, sampler, sobol_seed)
                               
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
              "'smooth!' will run Dykstra rounds to their floor on this domain where the matrix is dense.\n" *
              "Forcing pure alternating projections ('k_dykstra = 0') would need 'k_ap' of at least $(5*iters)."
    end
    if pred.rho_est > 0.99 && !pred.exact
        @warn "Parallel-plate-like coupling detected in the traced exchange factors.\n" *
              "ρ ≥ $(round(pred.rho_est; digits=6)) (lower bound, blind-set search truncated).\n" *
              "'smooth!' uses it for its Dykstra budget and densifies the matrix where it fits in memory;\n" *
              "only input too large to densify falls back to alternating projections, which may then be slow."
    end
    nothing
end