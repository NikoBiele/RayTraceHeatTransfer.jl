# ---------------- dual system: built once per mesh ----------------
struct DualSolver{T}
    Y::Matrix{T}        # reduced-mass weights (shared with OP)
    rowsum::Vector{T}   # Y*1  ⇒  R = Y + Diagonal(rowsum)
    dinv::Vector{T}     # Jacobi preconditioner 1 ./ diag(R)
end

function DualSolver(w::AbstractVector)
    Y = Y_mat(w)                                # threaded for N > 1000
    rowsum = vec(sum(Y, dims = 2))
    return DualSolver(Y, rowsum, 1 ./ (diag(Y) .+ rowsum))
end

Rmul!(y, S::DualSolver, v) = (mul!(y, Symmetric(S.Y, :U), v); @. y += S.rowsum * v; y)

function solve_R(S::DualSolver, b; rtol = 1e-14, maxiter = 200)
    x = zeros(eltype(b), length(b))
    bn = norm(b)
    bn == 0 && return x, 0                             # constraints already met exactly
    r = copy(b); z = S.dinv .* r; p = copy(z)
    rz = dot(r, z); Ap = similar(b)
    for it in 1:maxiter
        Rmul!(Ap, S, p)
        pAp = dot(p, Ap)
        pAp > 0 || return x, it                        # direction exhausted at rounding level
        alpha = rz / pAp
        @. x += alpha * p
        @. r -= alpha * Ap
        norm(r) ≤ rtol * bn && return x, it
        @. z = S.dinv * r
        rz_new = dot(r, z); beta = rz_new / rz
        @. p = z + beta * p
        rz = rz_new
    end
    @warn "solve_R: PCG reached maxiter = $maxiter, relative residual $(norm(r)/bn)"
    return x, maxiter
end

# convergence check
function delta_R(F::AbstractMatrix, w::Vector{T}) where T
    n = size(F, 1)
    d_squared = zero(T)
    @inbounds for i in 1:n
        w_i = w[i]
        for j in (i+1):n
            if F[i,j] > 1e-12 || F[j,i] > 1e-12
                w_j = w[j]
                x_ij = w_i * F[i,j]
                x_ji = F[j,i] * w_j
                diff = x_ij - x_ji
                d_squared += diff * diff / (w_i^2 + w_j^2)
            end
        end
    end
    return sqrt(d_squared)
end

# parallel convergence check
function delta_R_parallel(F::AbstractMatrix, w::Vector{T}) where T
    n = size(F, 1)
    d_squared = Threads.Atomic{T}(zero(T))
    Threads.@threads :dynamic for i in 1:n
        w_i = w[i]
        local_d_sum = zero(T)
        @inbounds for j in (i+1):n
            if F[i,j] > 1e-12 || F[j,i] > 1e-12
                w_j = w[j]
                x_ij = w_i * F[i,j]
                x_ji = F[j,i] * w_j
                diff = x_ij - x_ji
                local_d_sum += diff * diff / (w_i^2 + w_j^2)
            end
        end
        Threads.atomic_add!(d_squared, local_d_sum)
    end
    return sqrt(d_squared[])
end

function delta_R_X_sparse(X::SparseMatrixCSC, w::AbstractVector, u::Vector{Float64})
    N = size(X, 1)
    rows = rowvals(X); vals = nonzeros(X)
    nc = min(Threads.nthreads(), N)
    bounds = round.(Int, range(0, N, length = nc + 1))
    part = zeros(nc)
    @threads for c in 1:nc
        s = 0.0
        @inbounds for j in (bounds[c] + 1):bounds[c + 1]
            uj = u[j]; wj2 = w[j]^2
            @simd for k in nzrange(X, j)
                i = rows[k]
                if i < j
                    d = vals[k] * (u[i] - uj)
                    s += d * d / (w[i]^2 + wj2)
                end
            end
        end
        part[c] = s
    end
    return sqrt(sum(part))
end

function delta_R_X_dense(X::Matrix{Float64}, w::AbstractVector, u::Vector{Float64})
    N = size(X, 1)
    nc = min(4 * nthreads(), N)
    bounds = round.(Int, range(0, N, length = nc + 1))
    part = zeros(nc)
    @threads :dynamic for c in 1:nc
        s = 0.0
        @inbounds for j in (bounds[c] + 1):bounds[c + 1]
            uj = u[j]; wj2 = w[j]^2
            @simd for i in 1:(j - 1)
                d = X[i, j] * (u[i] - uj)
                s += d * d / (w[i]^2 + wj2)
            end
        end
        part[c] = s
    end
    return sqrt(sum(part))
end

# ---------------- raw reciprocity floor of the input (Lemma G.2 quantity) ----------------
@inline function _find_nz(F::SparseMatrixCSC, rows, i::Integer, j::Integer)
    r  = nzrange(F, j)
    lo = first(r); hi = last(r)
    @inbounds while lo <= hi
        m  = (lo + hi) >>> 1
        ri = rows[m]
        if ri == i
            return m
        elseif ri < i
            lo = m + 1
        else
            hi = m - 1
        end
    end
    return 0
end

function delta_R_raw(F::SparseMatrixCSC, w::AbstractVector)
    rows = rowvals(F)
    vals = nonzeros(F)
    n    = size(F, 2)
    T    = promote_type(eltype(vals), eltype(w))

    nt      = max(1, min(Threads.nthreads(), n))
    bounds  = round.(Int, range(1, n + 1; length = nt + 1))
    partial = zeros(T, nt)

    Threads.@threads for t in 1:nt
        s = zero(T)
        @inbounds for j in bounds[t]:(bounds[t+1] - 1)
            wj = w[j]
            for k in nzrange(F, j)
                i = rows[k]
                i == j && continue
                wi = w[i]
                if i < j
                    kt = _find_nz(F, rows, j, i)
                    d  = wi * vals[k] - (kt == 0 ? zero(T) : wj * vals[kt])
                else
                    _find_nz(F, rows, j, i) == 0 || continue
                    d = -(wi * vals[k])
                end
                s += d * d / (wi * wi + wj * wj)
            end
        end
        partial[t] = s
    end

    return sqrt(sum(partial))
end
delta_R_raw(F::AbstractMatrix, w::AbstractVector) =
    length(w) > 1000 ? delta_R_parallel(F, w) : delta_R(F, w)

function delta_perp(F::AbstractMatrix, w::AbstractVector;
                    mode_indicator::Union{Symbol,Nothing}=nothing,
                    dual::Union{Nothing,DualSolver}=nothing)
    N = length(w)
    if mode_indicator == :DYK                          # reciprocal iterate: δ_R = 0
        dual === nothing && (dual = DualSolver(w))
        b = w .* (vec(sum(F, dims = 2)) .- 1)
        lambda, _ = solve_R(dual, b)
        return sqrt(dot(b, lambda))
    end
    dR = delta_R_raw(F, w)
    mode_indicator == :AP && return dR                 # lower bound
    dual === nothing && (dual = DualSolver(w))
    _, b = Xbar_b(F, w, dual.Y)
    lambda, _ = solve_R(dual, b)
    return sqrt(dR^2 + dot(b, lambda))
end

function ap_distance_certificate(w::AbstractVector{<:Real};
                                 max_levels::Int = 64,
                                 bin_halfwidth::Real = 1e-3)
    N = length(w)
    a = sort!(unique(w))
    mu = if length(a) <= max_levels
        _mu_min_levels(a, _level_counts(w, a))               # exact
    elseif length(w) <= 2_000
        _mu_min_dense(w)                                     # exact, small N
    else                                                     # binned + margin
        ab, nb, dmax, dsum = _log_bin(w, bin_halfwidth)
        margin = 3 * (N * dmax + dsum) / 4                   # Weyl bound
        _mu_min_levels(ab, nb) - margin
    end
    mu = max(mu, 1.0)                                        # universal floor
    return sqrt(N / mu)
end

function _level_counts(w, a)
    n = zeros(Int, length(a))
    for wi in w
        n[searchsortedfirst(a, wi)] += 1
    end
    return n
end

function _mu_min_levels(a::AbstractVector{<:Real}, n::AbstractVector{<:Integer})
    L = length(a); N = sum(n)
    S = [a[l]^2 / (a[l]^2 + a[m]^2) for l in 1:L, m in 1:L]
    K = [a[l] * a[m] / (a[l]^2 + a[m]^2) for l in 1:L, m in 1:L]
    lam_within = S * n                       # eigenvalue per level, mult n_l - 1
    A = Symmetric(Diagonal(lam_within) - sqrt.(n * n') .* K)
    gmax = eigmax(A)
    for l in 1:L
        n[l] > 1 && (gmax = max(gmax, lam_within[l]))
    end
    return N - gmax
end

function _mu_min_dense(w::AbstractVector{<:Real})
    N = length(w)
    G = [-w[k] * w[l] / (w[k]^2 + w[l]^2) for k in 1:N, l in 1:N]
    for k in 1:N
        G[k, k] = sum(w[k]^2 / (w[k]^2 + w[j]^2) for j in 1:N if j != k)
    end
    return N - eigmax(Symmetric(G))
end

function _log_bin(w, t)
    lw = log.(w); lo = minimum(lw)
    nb = max(1, ceil(Int, (maximum(lw) - lo) / (2t)))
    bin = min.(nb, 1 .+ floor.(Int, (lw .- lo) ./ (2t)))
    s = zeros(nb); n = zeros(Int, nb)
    for (i, b) in enumerate(bin)
        s[b] += lw[i]; n[b] += 1
    end
    la = s ./ max.(n, 1)                       # bin means in log space
    d = abs.(lw .- la[bin])                    # per-element deviation
    keep = n .> 0
    return exp.(la[keep]), n[keep], maximum(d), sum(d)
end

function cross_coupling_chi(F::SparseMatrixCSC, n_surf::Integer; 
                            surfaces_first::Bool=true)
    N = size(F, 1)
    0 ≤ n_surf ≤ N  || throw(ArgumentError("n_surf = $n_surf out of range [0, $N]"))

    rows = rowvals(F)
    vals = nonzeros(F)
    acc  = zero(float(eltype(F)))

    # surface predicate for either ordering, lifted out of the hot path
    lo, hi = surfaces_first ? (1, n_surf) : (N - n_surf + 1, N)
    is_surf(idx) = lo ≤ idx ≤ hi

    @inbounds for j in 1:N
        col_surf = is_surf(j)
        for k in nzrange(F, j)
            i = rows[k]
            # keep the entry iff exactly one index is a surface ⇒ it lies in F_sg or F_gs
            (is_surf(i) ⊻ col_surf) && (acc += vals[k])
        end
    end
    chi = acc / N
    nz     = nnz(F)
    Ntot   = size(F, 1)
    perrow = nz / Ntot
    dens   = nz / Ntot^2
    mem_GB = Base.summarysize(F) / 2^30

    return chi, nz, perrow, dens, mem_GB
end

function cross_coupling_chi(F::AbstractMatrix, n_surf::Integer; surfaces_first::Bool=true)
    N = size(F, 1)
    0 ≤ n_surf ≤ N || throw(ArgumentError("n_surf = $n_surf out of range [0, $N]"))
    s = surfaces_first ? (1:n_surf) : ((N - n_surf + 1):N)
    g = surfaces_first ? ((n_surf + 1):N) : (1:(N - n_surf))
    acc = sum(@view F[s, g]) + sum(@view F[g, s])
    nz = count(!iszero, F)
    return acc / N, nz, nz / N, nz / N^2, Base.summarysize(F) / 2^30
end

# Gas escape fraction of a traced matrix: the fraction of gas emission absorbed
# at the walls, sum_g w_i sum_s F_ij / sum_g w_i, with w the emissive weights
# (areas for surfaces, 4 beta V for cells) and the surfaces first. Unlike the
# count-weighted chi, it is a property of the enclosure and converges under
# mesh refinement, so it serves as the coordinate of the AP / OP+AP rule.
function gas_escape_fraction(F::AbstractMatrix, w::AbstractVector, n_surf::Integer)
    N = size(F, 1)
    # No gas, or gas with no emissive weight (a transparent bin): the limit
    # beta -> 0, in which all gas emission reaches the walls.
    if n_surf >= N
        return 1.0
    end
    s = 1:n_surf
    g = (n_surf + 1):N
    w_g = @view w[g]
    total = sum(w_g)
    if total <= 0
        return 1.0
    end
    escape = vec(sum(F[g, s], dims = 2))
    return sum(w_g .* escape) / total
end

# Threshold of the gas escape fraction below which alternating projections
# alone give the more accurate solution, and above which the orthogonal
# projection start (OP+AP) does. Fitted on rectangular scattering enclosures
# against converged references, for 165 <= N <= 2805 elements and
# 10^3 to 10^4 rays per emitting element; both dependences are weak and
# logarithmic, and the fit is used as is outside that range.
function ap_threshold(N::Integer, rays_per_emitter::Real)
    return 0.583 - 0.081 * log10(N) - 0.024 * log10(rays_per_emitter)
end

function Y_mat(w::AbstractVector)
    N = length(w)
    if N > 1000
        Y_matrix = zeros(N,N)
        Threads.@threads for i in 1:N
            w2_i = w[i] * w[i]
            for j in 1:N
                w2_j = w[j] * w[j]
                Y_matrix[i,j] = (w2_i*w2_j)/(w2_i+w2_j)   # reduced-mass weights
            end
        end
        return Y_matrix
    else
        w2 = w .^ 2
        Y_matrix  = (w2 * w2') ./ (w2 .+ w2')                        # reduced-mass weights
        return Y_matrix
    end
end

# A priori bound on the Dykstra contraction rate from the support of the traced
# exchange factors. One Dykstra round spreads a row's residual over the row with
# the reduced-mass weights Y and the cone then discards the part that lands
# outside the support, so the one-row estimate of the contraction factor is
# 1 - (Y-weight on the support) / (total Y-weight). The symmetric coupling
# between rows can only slow this, so the slowest row gives a lower bound on
# the rate: the true rate is at least the value returned here.
function dykstra_rate_bound(F::AbstractMatrix, dual::DualSolver)
    N = size(F, 1)
    f_min = Inf
    for i in 1:N
        s = 0.0
        @inbounds for j in 1:N
            if F[i, j] > 0 || F[j, i] > 0
                s += dual.Y[i, j]
            end
        end
        f_min = min(f_min, s / dual.rowsum[i])
    end
    return 1 - f_min
end

function Xbar_b(F::AbstractMatrix, w::AbstractVector, Y::AbstractMatrix)
    N = length(w)
    if N > 1000
        Xbar = zeros(N, N); b = zeros(N); inv_w = 1 ./ w
        Threads.@threads for i in 1:N
            inv_wi = inv_w[i]; bi = 0.0
            @inbounds for j in 1:N
                v = Y[i, j] * (inv_wi * F[i, j] + F[j, i] * inv_w[j])
                Xbar[i, j] = v; bi += v
            end
            b[i] = bi - w[i]
        end
        return Xbar, b
    else
        Z = Diagonal(1 ./ w) * F
        Xbar = Y .* (Z + Z')
        return Xbar, Xbar * ones(N) - w
    end
end

function OP(F_OP::AbstractMatrix, w::AbstractVector, dual::DualSolver)
    Xbar, b = Xbar_b(F_OP, w, dual.Y)
    lambda, iters = solve_R(dual, b)
    Xstar = Xbar - dual.Y .* (lambda .+ lambda')      # rank-2 correction, one temporary
    return Diagonal(1 ./ w) * Xstar, iters
end

function DkAP(F_raw::AbstractMatrix, w::AbstractVector, num_surfaces::Int;
              k_dykstra::Int=0, k_ap::Int=1000, verbose::Bool=true,
              nz_over_N::Real=length(w), rho_ap::Union{Nothing,Real}=nothing)

    delta_raw = delta_R_raw(F_raw, w)
              
    if k_dykstra <= 0
        F_smooth, k_ap, converged, delta_max = AP(F_raw, w, num_surfaces; k_ap, verbose, nz_over_N)
        k_dykstra = 0; k_pcg_tot = 0; k_pcg_max = 0
        return F_smooth, converged, delta_raw, delta_max, k_dykstra, k_ap, k_pcg_tot, k_pcg_max
    end

    dual = DualSolver(w)                               # Y, R-diagonal, preconditioner: once
    N = length(w)
    F_smooth = copy(F_raw); P = zeros(size(F_raw)); delta = Inf; iters = 0
    k_pcg_tot = 0; k_pcg_max = 0

    # Cost rule. One Dykstra round is worth log(ρ_D)/log(ρ_AP) alternating
    # projection steps (the exchange rate) and costs (5 + k_pcg) N² against
    # 3 nnz(X) per step (the cost ratio). If a round pays, every round pays, so
    # Dykstra then runs to its floor; otherwise the prescribed k_dykstra rounds
    # are kept. ρ_AP comes from the a priori rate predictor, ρ_D from the
    # support bound; the measured δ⟂ ratios replace the bound as they arrive.
    k_cap = 20_000
    k_max = k_dykstra
    cost_ratio = Inf
    k_ap_needed = 0
    delta_hist = Tuple{Int,Float64}[]
    log_rho_ap = rho_ap === nothing ? nothing : log(min(float(rho_ap), 1.0))
    rho_d = dykstra_rate_bound(F_raw, dual)
    ap_cost = 3 * nz_over_N * N                         # 3·nnz of the traced matrix: the AP-only alternative
    target = 8 * eps(Float64)
    guard  = sqrt(N) * target                          # floor-aware acceptance, as in AP
    delta_prev = Inf; k_prev = 0
    delta_best = Inf; flat = 0
    k = 0
    while k < k_max
        k += 1
        G, iters = OP(F_smooth, w, dual)
        k_pcg_tot += iters
        k_pcg_max = max(iters, k_pcg_max)
        F_smooth = max.(G + P, 0.0)
        if k == 1 && log_rho_ap !== nothing
            cost_ratio = (5 + iters) * N^2 / ap_cost
            if log_rho_ap >= 0
                exchange = Inf
                k_ap_needed = typemax(Int)
            else
                exchange = log(rho_d) / log_rho_ap
                # alternating-projection steps from the raw defect to the floor at ρ_AP
                k_ap_needed = ceil(Int, log(guard / max(delta_raw, guard)) / log_rho_ap)
            end
            if exchange > cost_ratio
                k_max = k_cap
                verbose && println("Cost rule: exchange rate $(round(exchange; sigdigits=3)) > cost ratio $(round(cost_ratio; sigdigits=3)); Dykstra runs to its floor (ρ_AP ≈ $(round(rho_ap; sigdigits=6)), ρ_D ≥ $(round(rho_d; sigdigits=4)))")
            elseif k_ap_needed > k_ap
                k_max = k_cap
                verbose && println("Cost rule: exchange rate $(round(exchange; sigdigits=3)) ≤ cost ratio $(round(cost_ratio; sigdigits=3)), but alternating projections would need about $k_ap_needed steps (k_ap = $k_ap); Dykstra runs to its floor (ρ_AP ≈ $(round(rho_ap; sigdigits=6)), ρ_D ≥ $(round(rho_d; sigdigits=4)))")
            else
                verbose && println("Cost rule: exchange rate $(round(exchange; sigdigits=3)) ≤ cost ratio $(round(cost_ratio; sigdigits=3)); keeping $k_dykstra prescribed Dykstra rounds (ρ_AP ≈ $(round(rho_ap; sigdigits=6)), ρ_D ≥ $(round(rho_d; sigdigits=4)))")
            end
        end
        if k % 5 == 0 || k == k_max
            delta = delta_perp(F_smooth, w; mode_indicator=:DYK, dual)
            verbose && println("Dykstra round $k ($iters PCG iterations): δ⟂ = $delta")
            if delta >= delta_best * (1 - 1e-3)
                flat += 1
            else
                flat = 0
            end
            delta_best = min(delta_best, delta)
            if flat >= 3 && delta < guard
                verbose && println("  Dykstra contraction exhausted at round $k (floor δ⟂ ≈ $delta)")
                break
            end
            if isfinite(delta_prev) && delta_prev > 0 && delta > 0
                rho_meas = (delta / delta_prev)^(1 / (k - k_prev))
                # projected rounds to the floor, so a user can judge a long run early.
                # The local rate mostly rises as the slow mode takes over, so early
                # estimates tend to be low; the rate comes from a least-squares line
                # through log δ⟂ over the last four readings (20 rounds), printed
                # every 50 rounds. An estimate, not a bound.
                push!(delta_hist, (k, delta))
                if verbose && k_max == k_cap && k % 50 == 0 && delta > guard && length(delta_hist) >= 4
                    pts = delta_hist[end-3:end]
                    xs = [Float64(p[1]) for p in pts]
                    ys = [log(p[2]) for p in pts]
                    xm = sum(xs) / 4; ym = sum(ys) / 4
                    slope = sum((xs .- xm) .* (ys .- ym)) / sum((xs .- xm) .^ 2)
                    if slope < 0
                        rho_fit = exp(slope)
                        k_left = ceil(Int, log(guard / delta) / slope)
                        println("  measured ρ_D ≈ $(round(rho_fit; sigdigits=4)) over the last 20 rounds; estimated $k_left more rounds to the floor")
                    end
                end
                if log_rho_ap !== nothing && log_rho_ap < 0 && k_max == k_cap && k_ap_needed <= k_ap
                    if rho_meas >= 1
                        exchange_meas = 0.0
                    else
                        exchange_meas = log(rho_meas) / log_rho_ap
                    end
                    if exchange_meas < cost_ratio
                        verbose && println("Cost rule: measured exchange rate $(round(exchange_meas; sigdigits=3)) fell below the cost ratio; switching to alternating projections")
                        break
                    end
                end
            end
            delta_prev = delta; k_prev = k
        end
        # Below the guard sqrt(N)·8ε, continue to the rounding target 8ε only when
        # the measured rate makes that cheap (at most 100 more rounds); a slow run
        # stops at the guard, which is already sqrt(N) rounding units.
        if delta < target
            break
        end
        if delta < guard
            rounds_to_target = 100
            if isfinite(delta_prev) && delta_prev > delta && delta > 0 && k > k_prev
                rho_loc = (delta / delta_prev)^(1 / (k - k_prev))
                if rho_loc < 1
                    rounds_to_target = ceil(Int, log(target / delta) / log(rho_loc))
                end
            end
            if rounds_to_target >= 100
                break
            end
        end
        P = G + P - F_smooth
    end
    k_dyk_conv = k
    if k == k_cap && !(delta < max(target, guard) || (flat >= 3 && delta < guard))
        @warn "Dykstra reached the round cap $k_cap without reaching its floor; alternating projections finalise from here."
    end
    # A Dykstra iterate at its floor is already reciprocal and row-stochastic to
    # rounding, and δ⟂ is its exact distance to the feasible set. Returning it
    # directly keeps that certificate; the alternating-projection finaliser
    # cannot improve it and, near the parallel-plate limit, could only drift it
    # along the slow mode.
    F_smooth ./= sum(F_smooth, dims = 2)               # exact closure first; reciprocity stays at rounding level
    if delta < max(target, guard) || (flat >= 3 && delta < guard)
        verbose && println("Dykstra converged at round $k (δ⟂ = $delta); no finaliser needed")
        return F_smooth, true, delta_raw, delta, k_dyk_conv, 0, k_pcg_tot, k_pcg_max
    end
    # The finaliser accepts the same floor Dykstra used, sqrt(N)·8ε, so it does
    # not grind at the slow alternating-projection rate below a floor already reached.
    F_smooth, k_ap, converged, delta_max = AP(F_smooth, w, num_surfaces; k_ap, verbose, nz_over_N = 1.0)
    return F_smooth, converged, delta_raw, delta_max, k_dyk_conv, k_ap, k_pcg_tot, k_pcg_max
end

# DkAP restricted to the elements with positive weight. Elements with w == 0
# (volumes in a transparent bin, 4κV = 0) neither emit nor absorb in that bin:
# their columns are exactly zero and reciprocity with a zero weight is 0 = 0,
# so they are left out of the projection rather than regularised with a
# floor. Their rows, which the solver only ever multiplies by zero emission,
# are kept as traced and row-normalised over the active elements. Surfaces
# always have positive weight and come first, so num_surfaces is unchanged.
# Returns a matrix of the same size as F_raw.
function DkAP_active(F_raw::AbstractMatrix, w::AbstractVector, num_surfaces::Int;
                     k_dykstra::Int=0, k_ap::Int=1000, verbose::Bool=true,
                     nz_over_N::Real=length(w), rho_ap::Union{Nothing,Real}=nothing)
    active = findall(>(0), w)
    length(active) == length(w) &&
        return DkAP(F_raw, w, num_surfaces; k_dykstra, k_ap, verbose, nz_over_N, rho_ap)
    inactive = findall(<=(0), w)
    verbose && println("    $(length(inactive)) zero-weight elements excluded from the projection")

    F_act = F_raw[active, active]
    rs = sum(F_act, dims=2)
    rs[rs .== 0] .= 1
    F_act = F_act ./ rs
    out = DkAP(F_act, w[active], num_surfaces; k_dykstra, k_ap, verbose, nz_over_N, rho_ap)

    F_sub = out[1]
    F_smooth = F_sub isa SparseMatrixCSC ? spzeros(eltype(F_sub), size(F_raw)...) :
                                           zeros(eltype(F_sub), size(F_raw))
    F_smooth[active, active] = F_sub
    F_rows = F_raw[inactive, active]
    rs = sum(F_rows, dims=2)
    rs[rs .== 0] .= 1
    F_smooth[inactive, active] = F_rows ./ rs
    return (F_smooth, out[2:end]...)
end

function get_w(rtm::Union{RayTracingDomain2D,SurfaceDomain3D}; spectral_bin::Int=1)
    
    if typeof(rtm) <: RayTracingDomain2D
        # Extract surface areas, volumes, and local extinction for specific spectral bin
        w = zeros(length(rtm.surface_mapping)+length(rtm.volume_mapping))
        
        for ((coarse_idx, fine_idx, wall_idx), surface_idx) in rtm.surface_mapping
            face = rtm.fine_mesh[coarse_idx][fine_idx]
            w[surface_idx] = face.area[wall_idx]
        end
        for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
            face = rtm.fine_mesh[coarse_idx][fine_idx]
            element_idx = length(rtm.surface_mapping) + volume_idx
            # Extract local extinction for specific spectral bin
            if isa(face.kappa_g, Vector)
                local_beta = face.kappa_g[spectral_bin] + face.sigma_s_g[spectral_bin]
            else
                local_beta = face.kappa_g + face.sigma_s_g
            end
            w[element_idx] = 4*local_beta*face.volume
        end
        return w
    elseif typeof(rtm) <: SurfaceDomain3D
        w = Float64[]
        for emitter_face in rtm.facesMesh
            for emitter_sub in emitter_face.subFaces
                push!(w, emitter_sub.area)
            end
        end
        return w
    else
        error("Unknown domain type $rtm.")
    end
end

function get_b(rtm::Union{RayTracingDomain2D,SurfaceDomain3D})

    if typeof(rtm) <: RayTracingDomain2D
        b_mat = Matrix{Float64}(undef, length(rtm.surface_mapping)+length(rtm.volume_mapping), rtm.n_spectral_bins)    
        for m in 1:rtm.n_spectral_bins
            for ((coarse_idx, fine_idx, wall_idx), surface_idx) in rtm.surface_mapping
                face = rtm.fine_mesh[coarse_idx][fine_idx]
                b_mat[surface_idx,m] = 1 - face.epsilon[wall_idx][m]
            end
            for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
                face = rtm.fine_mesh[coarse_idx][fine_idx]
                element_idx = length(rtm.surface_mapping) + volume_idx
                β = face.sigma_s_g[m] + face.kappa_g[m]
                b_mat[element_idx,m] = β > 0 ? face.sigma_s_g[m] / β : 0.0   # transparent bin: no scattering
            end
        end
    elseif typeof(rtm) <: SurfaceDomain3D
        N_surfs = sum([length(superface.subFaces) for superface in rtm.facesMesh])
        b_mat = Matrix{Float64}(undef, N_surfs, rtm.n_spectral_bins)
        for m in 1:rtm.n_spectral_bins
            surf_count = 1
            for superface in rtm.facesMesh
                for face in superface.subFaces
                    b_mat[surf_count,m] = 1 - face.epsilon[m]
                    surf_count += 1
                end
            end
        end
    else
        error("Unknown domain type: $(typeof(rtm))")
    end
        
    return b_mat
end

function get_T_in(rtm::RayTracingDomain2D)
    # Extract surface areas, volumes, and local extinction for specific spectral bin
    T_in = Vector{Float64}(undef, length(rtm.surface_mapping)+length(rtm.volume_mapping))
    
    for ((coarse_idx, fine_idx, wall_idx), surface_idx) in rtm.surface_mapping
        face = rtm.fine_mesh[coarse_idx][fine_idx]
        T_in[surface_idx] = face.T_in_w[wall_idx]
    end
    for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
        face = rtm.fine_mesh[coarse_idx][fine_idx]
        element_idx = length(rtm.surface_mapping) + volume_idx
        T_in[element_idx] = face.T_in_g
    end
        
    return T_in
end

function get_T_current(rtm::RayTracingDomain2D)
    # Extract surface areas, volumes, and local extinction for specific spectral bin
    T_current = Vector{Float64}(undef, length(rtm.surface_mapping)+length(rtm.volume_mapping))
    
    for ((coarse_idx, fine_idx, wall_idx), surface_idx) in rtm.surface_mapping
        face = rtm.fine_mesh[coarse_idx][fine_idx]
        T_current[surface_idx] = face.T_w[wall_idx]
    end
    for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
        face = rtm.fine_mesh[coarse_idx][fine_idx]
        element_idx = length(rtm.surface_mapping) + volume_idx
        T_current[element_idx] = face.T_g
    end
        
    return T_current
end

# Updated smoothing algorithm
function smooth_F(rtm::Union{RayTracingDomain2D,SurfaceDomain3D}, F_raw::AbstractMatrix;
                                 k_ap::Int=1_000,
                                 smooth_surfaces_only::Bool=false,
                                 k_dykstra::Int=1,
                                 verbose::Bool=true,
                                 renorm::Bool=true,
                                 spectral_bin::Int=1,
                                 rho_ap::Union{Nothing,Real}=nothing)

    w = get_w(rtm; spectral_bin)
    num_surfaces = if typeof(rtm) <: RayTracingDomain2D
        length(rtm.surface_mapping)
    else
        size(F_raw,1)
    end
    N = length(w)

    verbose && println("Matrix size: $(N)×$(N)")

    # Nonzero statistics decide the sparse/dense path; the gas escape fraction
    # Φ_g against its threshold decides AP-only versus OP+AP (2D participating
    # media only). The threshold depends weakly on N and on the rays per
    # emitting element; a domain that was never traced reports n_rays = 0
    # and gets the 10^4 rays per emitter the rule was centred on.
    nz = N^2 # initialize dense (for view factors)
    phi_g = 1.0      # no gas: the beta -> 0 limit, in which all gas emission reaches the walls
    phi_star = 0.0
    if rtm.surfaces_only
        if F_raw isa SparseMatrixCSC
            nz = nnz(F_raw)
        end
    else
        chi, nz, perrow, dens, mem_GB = cross_coupling_chi(F_raw, num_surfaces)
        verbose && println("Raw ray traced sparse matrix statistics:")
        verbose && println("    nonzeros: ", nz, ", per-row density: ", round(perrow; sigdigits=6),
                            ", density: ", round(dens;sigdigits=6), ", memory footprint (GB): ", mem_GB)
        phi_g = gas_escape_fraction(F_raw, w, num_surfaces)
        if !isfinite(phi_g)
            phi_g = 1.0
        end
        if rtm.n_rays > 0
            rays_per_emitter = rtm.n_rays / N
        else
            rays_per_emitter = 10_000.0
        end
        phi_star = ap_threshold(N, rays_per_emitter)
        verbose && println("Gas escape fraction of F_raw is Φ_g = $(round(phi_g; sigdigits=4)); AP threshold Φ_g* = $(round(phi_star; sigdigits=4)) at N = $N, $(round(Int, rays_per_emitter)) rays per emitter")
    end
    nz_over_N = Float64(nz / N)
    if nz / N^2 > 0.25
        F_raw = Matrix(F_raw)
    end
    ap_preferred = phi_g < phi_star

    # Near-parallel-plate input (ρ_AP > 0.99) overrides the Φ_g rule, and the
    # sparse gate when the dense Dykstra arrays fit in memory, because pure
    # alternating projections are then infeasible; the cost rule in DkAP sets
    # the rounds. Dykstra holds about six dense N×N Float64 arrays.
    # Sparse input above the threshold is densified for Dykstra when the dense
    # working set (about eight N×N Float64 arrays, 64 N² bytes) stays within
    # ~1 GB, i.e. N up to about 4000; near-parallel-plate input is densified
    # whenever the arrays fit in memory, because pure alternating projections
    # are infeasible there.
    flagged     = rho_ap !== nothing && rho_ap > 0.99
    dense_fits  = 64 * N^2 < 0.5 * Sys.free_memory()
    dense_small = 64 * N^2 < 2^30
    densify     = F_raw isa SparseMatrixCSC && ((flagged && dense_fits) || (!ap_preferred && dense_small))
    if typeof(rtm) <: ViewFactorDomain3D ||
        (typeof(rtm) <: RayTracingDomain2D && ap_preferred && !flagged) ||
        (F_raw isa SparseMatrixCSC && !densify)

        if flagged && F_raw isa SparseMatrixCSC
            verbose && println("    Near-parallel-plate input too large to densify; using AP only")
        else
            verbose && println("    Using AP only (0 Dykstra rounds)")
        end
        k_dykstra = 0

    else
        if F_raw isa SparseMatrixCSC
            verbose && println("    Sparse input above the AP threshold with N = $N: densifying for Dykstra")
            F_raw = Matrix(F_raw)
        end
        if rho_ap === nothing
            verbose && println("    Using prescribed $k_dykstra Dykstra rounds")
        else
            verbose && println("    Dykstra budget set by the cost rule (prescribed k_dykstra = $k_dykstra applies if a round does not pay)")
        end
    end

    if smooth_surfaces_only && !(F_raw isa SparseMatrixCSC)
        w = renorm ? w[1:num_surfaces]./minimum(w[1:num_surfaces]) : w[1:num_surfaces]
        F_raw = F_raw[1:num_surfaces,1:num_surfaces]
    else
        w = renorm ? w./minimum(w[w .> 0]) : w
    end

    return DkAP_active(F_raw, w, num_surfaces; k_dykstra=k_dykstra,
                k_ap=k_ap, verbose=verbose, nz_over_N=nz_over_N, rho_ap=rho_ap)
end

function AP_convergence_check(w::AbstractVector, num_surfaces::Int)
    num_volumes = num_surfaces == length(w) ? 0 : length(w) - num_surfaces
    if num_volumes == 0
        w_max = maximum(w)
        w_sum = sum(w)
        check_passed = w_max < 0.5 * w_sum
        if !check_passed
            error("Smoothing convergence check failed: max surface w ($w_max) ≥ half of total w ($(w_sum/2));\n
                   The enclosure is physically impossible and smoothing convergence not possible.")
        end
    end
end

function build_X(F::SparseMatrixCSC, w::AbstractVector)
    X = Diagonal(w) * F                  # sparse, scales rows
    return 0.5*(X + X')  # sparse + sparse, O(nnz)
end

function build_X(F::Matrix{Float64}, w::Vector{Float64})
    N = size(F, 1)
    X = similar(F)
    @threads for j in 1:N
        wj = w[j]
        @inbounds @simd for i in 1:N
            X[i, j] = 0.5 * (w[i] * F[i, j] + wj * F[j, i])
        end
    end
    return X
end

# row sums and hunger vector of the current iterate
function hunger!(r::Vector{Float64}, u::Vector{Float64}, X::Matrix{Float64},
                 w::AbstractVector, e::Vector{Float64})
    mul!(r, Symmetric(X, :U), e)
    u .= w ./ r
end

function hunger!(r::Vector{Float64}, u::Vector{Float64}, X::SparseMatrixCSC,
                 w::AbstractVector, ::Vector{Float64})
    vals = nonzeros(X)
    @threads for j in 1:size(X, 1)                     # column sums = row sums
        s = 0.0
        @inbounds @simd for k in nzrange(X, j)
            s += vals[k]
        end
        r[j] = s; u[j] = w[j] / s
    end
    return u
end

# Hadamard scaling X_k -> X_{k+1}
function scale!(X::Matrix{Float64}, u::Vector{Float64})
    @threads for j in 1:size(X, 1)
        uj = u[j]
        @inbounds @simd for i in 1:size(X, 1)
            X[i, j] *= 0.5 * (u[i] + uj)
        end
    end
    return X
end

function scale!(X::SparseMatrixCSC, u::Vector{Float64})
    rows = rowvals(X); vals = nonzeros(X)
    @threads for j in 1:size(X, 1)
        uj = u[j]
        @inbounds @simd for k in nzrange(X, j)
            vals[k] *= 0.5 * (u[rows[k]] + uj)
        end
    end
    return X
end

delta_R_X(X::SparseMatrixCSC, w, u) = delta_R_X_sparse(X, w, u)
delta_R_X(X::Matrix{Float64}, w, u) = delta_R_X_dense(X, w, u)

# r is the current row-sum vector of X (from the last hunger! call); no recomputation
function recover_F(X::SparseMatrixCSC, r::Vector{Float64})
    rows = rowvals(X); vals = nonzeros(X)
    @threads for j in 1:size(X, 1)
        @inbounds @simd for k in nzrange(X, j)
            vals[k] /= r[rows[k]]
        end
    end
    return X
end
recover_F(X::Matrix{Float64}, r::Vector{Float64}) = X ./ r

function AP(F_AP::AbstractMatrix, w::AbstractVector, num_surfaces::Int;
            k_ap::Int=1_000, verbose::Bool=true,
            nz_over_N::Real=length(w))

    AP_convergence_check(w, num_surfaces)
    mult = ap_distance_certificate(w)
    use_sparse = F_AP isa SparseMatrixCSC
    use_sparse || (F_AP = Matrix(F_AP))
    verbose && println("Alternating projection (AP): parallel, $(use_sparse ? "sparse" : "dense")")

    N = length(w)
    target = 8 * eps(Float64)
    guard  = sqrt(N / nz_over_N) * target              # floor-aware acceptance: 8·N/√nnz·ε, = 8ε when dense
    switch = 4 * guard
    max_stride = 52

    delta_raw = delta_R_raw(F_AP, w)                   # defect of the input (MC-bounded)
    X_iter = build_X(F_AP, w)
    r = zeros(N); u = zeros(N); e = ones(N)
    hunger!(r, u, X_iter, w, e)                        # (X_0, u_0)
    delta = delta_R_X(X_iter, w, u)                    # defect of the iterate recovered from X_0
    delta < target && (k_ap = 0)                    # already converged, exit
    delta_init = delta; delta_best = delta
    verbose && println("  raw input: δ_R = $delta_raw;  after first reciprocity projection: $delta ≤ δ⟂ ≤ $(mult*delta)")

    k = 0; k_next = 0; c = 0; flat = 0
    k_prev = 0; delta_prev = delta; rho_est = 0.5; floor_accepted = false

    while k < k_ap && delta > target
        scale!(X_iter, u)                              # X_k -> X_{k+1}
        k += 1
        hunger!(r, u, X_iter, w, e)                    # (X_k, u_k), needed for the next step anyway
        if k >= k_next
            delta = delta_R_X(X_iter, w, u)            # δ_R of the current iterate, no lag
            c += 1
            if c >= 3
                rho_est = min(max((delta / delta_prev)^(1 / max(k - k_prev, 1)), 0.5), 0.9999)
                flat = delta >= delta_best * (1 - 1e-3) ? flat + 1 : 0
            end
            delta_best = min(delta_best, delta)
            if delta < guard && (rho_est > 0.99 || flat >= 3)
                floor_accepted = true
                verbose && println("  contraction exhausted at iteration $k (floor δ_R ≈ $delta)")
                break
            end
            k_prev, delta_prev = k, delta
            if delta > switch
                stride = max(1, ceil(Int, log(delta / target) / log(1 / rho_est)))
                k_next = k + min(stride, max_stride)
            else
                k_next = k + 1
            end
            verbose && println("  iteration $k: $delta ≤ δ⟂ ≤ $(mult*delta)")
        end
    end
    converged = delta <= target || floor_accepted
    if converged
        verbose && println("Converged after $k iterations. Final: $delta ≤ δ⟂ ≤ $(mult*delta)")
    else
        @warn "AP reached k_ap = $k_ap. Final: $delta ≤ δ⟂ ≤ $(mult*delta)"
    end
    if delta > max(delta_init, guard)
        @warn "Smoothing increased the distance to the target manifold; use F_raw instead of F_smooth."
    end
    return recover_F(X_iter, r), k, converged, (mult*delta)
end