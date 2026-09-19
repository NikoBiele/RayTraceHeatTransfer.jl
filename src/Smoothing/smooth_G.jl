# Smoothing of the angular exchange factors G[a][i, j].
#
# Constraints: G ≥ 0 with its zero pattern kept; per-bin row sums
# Σ_j G[a][i, j] = p_i[a] (the analytic emission shares, which the tracer samples
# by construction); reversed-bin reciprocity w_i G[a][i, j] = w_j G[ā][j, i].
# On the (element, bin) nodes this is the ordinary problem for the row-stochastic
# matrix 𝐅[(i,a),(j,ā)] = G[a][i, j] / p_i[a] with weights w_i p_i[a], and it
# splits into independent problems, one per reversal pair {a, ā}: a bipartite
# 2N × 2N block [[0, F_a], [F_ā, 0]]. Only the Hadamard-scaling alternating
# projection is used: it preserves the zero pattern, which the Dykstra/OP
# correction (a dense rank-2 update) would fill in.

# one reversal pair; nodes 1:N ↔ (·, a), N+1:2N ↔ (·, ā)
function _smooth_G_pair(Ga::SparseMatrixCSC, Gā::SparseMatrixCSC, w::AbstractVector,
                        pa::AbstractVector, pā::AbstractVector, ns::Int; k_ap::Int)
    N = size(Ga, 1)
    Is = Int[]; Js = Int[]; Vs = Float64[]
    Ia, Ja, Va = findnz(Ga)
    for t in eachindex(Va)
        i = Ia[t]
        (pa[i] > 0 && Va[t] > 0) || continue
        push!(Is, i); push!(Js, N + Ja[t]); push!(Vs, Va[t] / pa[i])
    end
    Ib, Jb, Vb = findnz(Gā)
    for t in eachindex(Vb)
        i = Ib[t]
        (pā[i] > 0 && Vb[t] > 0) || continue
        push!(Is, N + i); push!(Js, Jb[t]); push!(Vs, Vb[t] / pā[i])
    end
    Fp = sparse(Is, Js, Vs, 2N, 2N)
    wp = vcat(w .* pa, w .* pā)
    if !any(>(0), wp) || nnz(Fp) == 0                      # nothing to project (e.g. fully transparent)
        return spzeros(N, N), spzeros(N, N), true, 0.0, 0.0, 0
    end
    wp ./= minimum(wp[wp .> 0])
    Fs, converged, delta_raw, delta_max, _, k_used, _, _ =
        DkAP_active(Fp, wp, ns; k_dykstra = 0, k_ap = k_ap, verbose = false,
                    nz_over_N = nnz(Fp) / (2N))
    Gsa = sparse(Diagonal(pa) * Fs[1:N, N+1:2N])
    Gsā = sparse(Diagonal(pā) * Fs[N+1:2N, 1:N])
    return Gsa, Gsā, converged, delta_raw, delta_max, k_used
end

function smooth_G!(rtm::RayTracingDomain2D; k_ap::Int = 20_000, verbose::Bool = true, keep_F_raw::Bool = true)
    dm = rtm.directional_model
    dm === nothing && error("smooth_G!: the domain has no directional_model")
    G_raw = rtm.G_raw
    G_raw === nothing && error("no angular exchange factors: set directional_model, then trace with " *
                               "domain(N_rays; method = :pathlength)")
    spectral_layout = G_raw isa Vector{<:Vector}
    K  = spectral_layout ? length(G_raw) : 1
    A  = n_angular_bins(dm)
    ns = length(rtm.surface_mapping)
    p  = emission_share_matrix(rtm)
    pair_list = [a for a in 1:A if a < reversed_bin(dm, a)]

    G_all = Vector{Vector{SparseMatrixCSC{Float64,Int}}}(undef, K)
    F_all = Vector{AbstractMatrix}(undef, K)
    converged_out = Vector{Bool}(undef, K);    delta_raw_out = Vector{Float64}(undef, K)
    delta_max_out = Vector{Float64}(undef, K); k_ap_out      = Vector{Int}(undef, K)

    for k in 1:K
        Gk = spectral_layout ? G_raw[k] : G_raw
        w  = get_w(rtm; spectral_bin = k)
        Gs = Vector{SparseMatrixCSC{Float64,Int}}(undef, A)
        conv = trues(length(pair_list)); d_raw = zeros(length(pair_list))
        d_max = zeros(length(pair_list)); k_it = zeros(Int, length(pair_list))
        Threads.@threads for n in eachindex(pair_list)
            a = pair_list[n]
            ā = reversed_bin(dm, a)
            Gs[a], Gs[ā], conv[n], d_raw[n], d_max[n], k_it[n] =
                _smooth_G_pair(Gk[a], Gk[ā], w, view(p, :, a), view(p, :, ā), ns; k_ap = k_ap)
        end
        G_all[k] = Gs
        Fk = reduce(+, Gs)
        F_all[k] = nnz(Fk) / length(Fk) > 0.25 ? Matrix(Fk) : Fk   # same rule as smooth_F
        converged_out[k] = all(conv);        delta_raw_out[k] = maximum(d_raw)
        delta_max_out[k] = maximum(d_max);   k_ap_out[k]      = maximum(k_it)
        verbose && println("Smoothed G for spectral bin $k/$K: $(length(pair_list)) reversal pairs, ",
                           "AP iterations ≤ $(k_ap_out[k]), δ_R raw ≤ $(round(delta_raw_out[k], sigdigits = 3)), ",
                           "bound after ≤ $(round(delta_max_out[k], sigdigits = 3))",
                           converged_out[k] ? "" : "  — NOT converged, raise k_ap")
    end

    rtm.G_smooth = spectral_layout ? G_all : G_all[1]
    rtm.F_smooth = rtm.F_raw isa Vector ? F_all : F_all[1]   # F_smooth = Σₐ G_smooth[a]
    keep_F_raw || (rtm.F_raw = Matrix{Float64}(undef, 2, 2))
    return (; converged = converged_out, delta_raw = delta_raw_out, delta_smooth = delta_max_out,
              k_dykstra = zeros(Int, K), k_ap = k_ap_out, k_pcg_tot = zeros(Int, K), k_pcg_max = zeros(Int, K))
end