struct DualSolver{T}
    M::Matrix{T}        # reduced-mass weights (OP needs these anyway)
    rowsum::Vector{T}   # M*1  ⇒  R = M + Diagonal(rowsum)
    dinv::Vector{T}     # 1 ./ diag(R), Jacobi preconditioner
end

function DualSolver(w::AbstractVector)
    w2 = w .^ 2
    M = (w2 * w2') ./ (w2 .+ w2')
    rowsum = vec(sum(M, dims = 2))
    return DualSolver(M, rowsum, 1 ./ (diag(M) .+ rowsum))
end

Rmul!(y, S::DualSolver, v) = (mul!(y, S.M, v); @. y += S.rowsum * v; y)

function solve_R(S::DualSolver, b; rtol = 1e-14, maxiter = 200)
    x = zeros(eltype(b), length(b))
    r = copy(b); z = S.dinv .* r; p = copy(z)
    rz = dot(r, z); bn = norm(b); Ap = similar(b)
    for it in 1:maxiter
        Rmul!(Ap, S, p)
        alpha = rz / dot(p, Ap)
        @. x += alpha * p
        @. r -= alpha * Ap
        norm(r) ≤ rtol * bn && return x, it
        @. z = S.dinv * r
        rz_new = dot(r, z); beta = rz_new / rz
        @. p = z + beta * p
        rz = rz_new
    end
    return x, maxiter
end

# steps 1-4: F = 1/2*( F + E^(-1)*F^T*E )/sum( 1/2 + 1/2*E^(-1)*F^T*E )
function all_steps!(F::AbstractMatrix, E::Vector{T}) where T
    n = size(F, 1)
    @inbounds for i in 1:n
        E_i = E[i]
        inv_E_i = 1/E_i
        for j in (i+1):n
            if F[i,j] > 1e-12 || F[j,i] > 1e-12
                F_ij = F[i,j]
                F_ji = F[j,i]
                E_j = E[j]
                EF = 0.5 * (E_i * F_ij + F_ji * E_j)
                F[i,j] = inv_E_i * EF
                F[j,i] = (1 / E_j) * EF
            end
        end
    end
    @inbounds for i in 1:n
        row_sum = zero(T)
        for j in 1:n
            row_sum += F[i,j] # Compute row sum
        end
        # Normalize row
        inv_sum = 1.0 / row_sum
        for j in 1:n
            F[i, j] *= inv_sum
        end
    end
end

# convergence check
function delta_R(F::AbstractMatrix, E::Vector{T}) where T
    n = size(F, 1)
    d_squared = zero(T)
    @inbounds for i in 1:n
        E_i = E[i]
        for j in (i+1):n
            if F[i,j] > 1e-12 || F[j,i] > 1e-12
                F_ij = F[i,j]
                F_ji = F[j,i]
                E_j = E[j]
                x_ij = E_i * F_ij
                x_ji = F_ji * E_j
                diff = x_ij - x_ji
                d_squared += diff * diff / (E_i^2 + E_j^2)
            end
        end
    end
    return sqrt(d_squared)
end

# steps 1-4: F = 1/2*( F + E^(-1)*F^T*E )/sum( 1/2 + 1/2*E^(-1)*F^T*E )
function all_steps_parallel!(F::AbstractMatrix, E::Vector{T}) where T
    n = size(F, 1)
    Threads.@threads :dynamic for i in 1:n
        E_i = E[i]
        inv_E_i = 1/E_i
        @inbounds for j in (i+1):n
            if F[i,j] > 1e-12 || F[j,i] > 1e-12
                F_ij = F[i,j]
                F_ji = F[j,i]
                E_j = E[j]
                EF = 0.5 * (E_i * F_ij + F_ji * E_j)
                F[i,j] = inv_E_i * EF
                F[j,i] = (1 / E_j) * EF
            end
        end
    end
    Threads.@threads for i in 1:n
        row_sum = zero(T)
        @inbounds for j in 1:n
            row_sum += F[i,j] # Compute row sum
        end
        # Normalize row
        inv_sum = 1.0 / row_sum
        @inbounds for j in 1:n
            F[i, j] *= inv_sum
        end
    end
end

# parallel convergence check
function delta_R_parallel(F::AbstractMatrix, E::Vector{T}) where T
    n = size(F, 1)
    d_squared = Threads.Atomic{T}(zero(T))
    Threads.@threads :dynamic for i in 1:n
        E_i = E[i]
        local_d_sum = zero(T)
        @inbounds for j in (i+1):n
            if F[i,j] > 1e-12 || F[j,i] > 1e-12
                F_ij = F[i,j]
                F_ji = F[j,i]
                E_j = E[j]
                x_ij = E_i * F_ij
                x_ji = F_ji * E_j
                diff = x_ij - x_ji
                local_d_sum += diff * diff / (E_i^2 + E_j^2)
            end
        end
        Threads.atomic_add!(d_squared, local_d_sum)
    end
    return sqrt(d_squared[])
end

function delta_perp(F::AbstractMatrix, w::AbstractVector;
                    mode_indicator::Union{Symbol,Nothing}=nothing,
                    tperm::Union{Nothing,Vector{Int}}=nothing)
    N = length(w)
    if mode_indicator == :DYK
        # Dykstra rounds (reciprocity violation = 0)
        b = Diagonal(w)*(F*ones(N)-ones(N))
        # lambda = R_chol \ b
        dual = DualSolver(w)
        lambda, iters = solve_R(dual, b)
        return sqrt(b' * lambda)
    else
        if F isa SparseMatrixCSC
            delta_R_val = reciprocity_residual_sparse(F, w, tperm)
        elseif N > 1000
            delta_R_val = delta_R_parallel(F, w)
        else
            delta_R_val = delta_R(F, w)
        end
        if mode_indicator == :AP_sparse
            return delta_R_val # estimate (lower bound)
        else # mode_indicator == :AP_dense
            # full distance
            dual = DualSolver(w)
            _ , b = Xbar_b(F, w, dual.M)
            lambda, iters = solve_R(dual, b)
            return sqrt(delta_R_val^2 + b' * lambda)
        end
    end
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
    gmax = maximum(eigvals(A))
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
    return N - maximum(eigvals(Symmetric(G)))
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

# sparse matrices

function symmetric_workspace(F::SparseMatrixCSC{Tv}) where {Tv}
    n = size(F, 2)
    W = F + F' + sparse(I, n, n)         # symmetric pattern + stored diagonal (for build_M!)
    fill!(nonzeros(W), zero(Tv))
    rows = rowvals(F); vals = nonzeros(F)
    @inbounds for j in 1:n, idx in nzrange(F, j)
        W[rows[idx], j] = vals[idx]
    end
    return W, build_tperm(W)
end

# tperm[k] = storage slot of the transpose of slot k. Requires structural symmetry.
function build_tperm(F::SparseMatrixCSC)
    n = size(F,2); size(F,1) == n || error("square required")
    rows = rowvals(F)
    cp = SparseArrays.getcolptr(F)
    cursor = collect(cp[1:n])         # next unconsumed slot per column
    tperm = Vector{Int}(undef, nnz(F))
    for j in 1:n, idx in nzrange(F, j)
        i = rows[idx]
        t = cursor[i]
        rows[t] == j || error("not structurally symmetric at ($i,$j)")
        tperm[idx] = t
        cursor[i] += 1
    end
    return tperm
end

function all_steps_resid_sparse!(F::SparseMatrixCSC, E::AbstractVector, tperm::Vector{Int},
                          tmp::Vector, rowsum::Vector)
    all_steps_sparse!(F, E, tperm, tmp, rowsum)
    return reciprocity_residual_sparse(F, E, tperm)
end

# ‖ Diagonal(E)·F − Fᵀ·Diagonal(E) ‖_F   (reciprocity defect Eᵢ Fᵢⱼ − Eⱼ Fⱼᵢ)
function reciprocity_residual_sparse(F::SparseMatrixCSC, E::AbstractVector, tperm::Vector{Int})
    rows = rowvals(F); vals = nonzeros(F)
    acc = Threads.Atomic{eltype(vals)}(zero(eltype(vals)))
    Threads.@threads for k in eachindex(vals)
        E_k  = E[rows[k]]
        E_kT = E[rows[tperm[k]]]
        r = E_k * vals[k] - E_kT * vals[tperm[k]]
        Threads.atomic_add!(acc, r * r / (E_k^2 + E_kT^2))
    end
    return sqrt(acc[] / 2)
end

function cross_coupling_chi(F::SparseMatrixCSC, n_surf::Integer; surfaces_first::Bool=true)
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

function all_steps_sparse!(F::SparseMatrixCSC, E::AbstractVector, tperm::Vector{Int},
                    tmp::Vector, rowsum::Vector)
    rows = rowvals(F); vals = nonzeros(F)

    Threads.@threads for k in eachindex(vals)          # 1) Diagonal(E) * F
        vals[k] *= E[rows[k]]
    end

    copyto!(tmp, vals)                          # 2) 0.5*(F + F')
    Threads.@threads for k in eachindex(vals)
        vals[k] = 0.5 * (tmp[k] + tmp[tperm[k]])
    end

    Threads.@threads for k in eachindex(vals)          # 3) Diagonal(1/E) * F
        vals[k] /= E[rows[k]]
    end

    fill!(rowsum, zero(eltype(rowsum)))         # 4) row-normalize
    @inbounds for k in eachindex(vals)
        rowsum[rows[k]] += vals[k]
    end
    Threads.@threads for k in eachindex(vals)
        vals[k] /= rowsum[rows[k]]
    end
    return F
end

function M_mat(w::AbstractVector)
    N = length(w)
    if N > 1000
        M = zeros(N,N)
        Threads.@threads for i in 1:N
            w2_i = w[i] * w[i]
            for j in 1:N
                w2_j = w[j] * w[j]
                M[i,j] = (w2_i*w2_j)/(w2_i+w2_j)   # reduced-mass weights
            end
        end
        return M
    else
        w2 = w .^ 2
        M  = (w2 * w2') ./ (w2 .+ w2')                        # reduced-mass weights
        return M
    end
end

function Xbar_b(F::AbstractMatrix, w::AbstractVector, M::AbstractMatrix)
    N = length(w)
    if N > 1000
        Xbar = zeros(N,N)
        b = zeros(N)
        Threads.@threads for i in 1:N
            inv_wi = 1 / w[i]
            for j in 1:N
                Xbar[i,j] = M[i,j] * (inv_wi * F[i,j] + F[j,i] * (1 ./ w[j]))  # weighted symmetrisation
                b[i] += Xbar[i,j]
            end
            b[i] -= w[i] # conservation defect
        end
        return Xbar, b
    else
        Z  = Diagonal(1 ./ w) * F                            # W⁻² X0
        Xbar  = M .* (Z + Z')                                     # weighted symmetrisation
        b  = Xbar * ones(N) - w                                   # conservation defect
        return Xbar, b
    end
end

function OP(F::AbstractMatrix, w::AbstractVector;
            M::Union{Nothing,AbstractMatrix}=nothing) # where P
    N = length(w)
    Xbar, b = Xbar_b(F, w, M)
    dual = DualSolver(w)
    lambda, iters = solve_R(dual, b)
    Xstar = Xbar - M .* (lambda * ones(N)' + ones(N) * lambda')            # rank-2 correction
    return Diagonal(1 ./ w) * Xstar, iters
end

function DkAP(F_raw::AbstractMatrix, w::AbstractVector;
              k_dykstra::Union{Nothing,Int}=nothing, 
              max_iters::Int=1000, verbose::Bool=true,
              nz_over_N::AbstractFloat=1.0)

    if k_dykstra > 0
        F_smooth = copy(F_raw); P = zeros(size(F_raw))
        M = M_mat(w); delta = 100.0; k = 0
        for k in 1:k_dykstra
            G, iters = OP(F_smooth, w; M=M)
            F_smooth = max.(G + P, 0.0)
            if k % 5 == 0
                delta = delta_perp(F_smooth, w; mode_indicator=:DYK)
            end
            if delta < 8*eps(Float64)
                verbose && println("Dykstra converged in $k rounds ($iters iterations in final inner solve)")
                return F_smooth
            end
            if k < k_dykstra
                P = G + P - F_smooth
            end
        end
        F_smooth = F_smooth ./ sum(F_smooth, dims=2)
        return AP(F_smooth, w; max_iters=max_iters,
                    nz_over_N=nz_over_N,
                    verbose=verbose, M=M)
    else
        mult = ap_distance_certificate(w)
        return AP(F_raw, w; max_iters=max_iters, verbose=verbose,
                    nz_over_N=nz_over_N, mult=mult)
    end
end

function get_w(rtm::RayTracingDomain2D; spectral_bin::Int=1)
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
        w[element_idx] = max(1e-6, 4*local_beta*face.volume)
    end
    
    return w
end

function get_b(rtm::Union{RayTracingDomain2D,ViewFactorDomain3D})

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
                b_mat[element_idx,m] = face.sigma_s_g[m]/(face.sigma_s_g[m]+face.kappa_g[m])
            end
        end
    elseif typeof(rtm) <: ViewFactorDomain3D
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
function smooth_F(F_raw::AbstractMatrix, w::AbstractVector;
                                 max_iters::Int=1_000,
                                 smooth_surfaces_only::Bool=false,
                                 k_dykstra::Union{Nothing,Int}=nothing,
                                 verbose::Bool=true,
                                 renorm::Bool=true)

    # count zero diagonal (zero self-view for planar surfaces)
    num_surfaces = sum([F_raw[i,i] < 4*eps(Float64) for i in 1:length(w)])

    # use for deciding between sparse/dense path and AP/OP
    if smooth_surfaces_only
        off_diag_chi = 0.0 # convex enclosure (deterministic exchange factor, only need polishing, i.e. AP)
        nz_over_N = Float64(length(w))
    else
        off_diag_chi, nz, perrow, dens, mem_GB = cross_coupling_chi(F_raw, num_surfaces)
        verbose && println("Raw ray traced sparse matrix statistics:")
        verbose && println("    nonzeros: ", nz, ", per-row density: ", round(perrow; sigdigits=6),
                            ", density: ", round(dens;sigdigits=6), ", memory footprint (GB): ", mem_GB)
        nz_over_N = Float64(nz / length(w))
        if nz / length(w)^2 > 0.25
            F_raw = Matrix(F_raw)
        end
    end

    verbose && println("Off-diagonal cross-coupling of F_raw is χ = $off_diag_chi")
    if k_dykstra === nothing
        if off_diag_chi < 0.4 || F_raw isa SparseMatrixCSC
            verbose && println("    Using AP only (0 Dykstra rounds)")
            k_dykstra = 0
        else
            verbose && println("    Using OP+AP (1 Dykstra round)")
            k_dykstra = 1
        end
    else
        verbose && println("    Using prescribed $k_dykstra Dykstra rounds")
    end

    if smooth_surfaces_only && !(F_raw isa SparseMatrixCSC)
        w = renorm ? w[1:num_surfaces]./minimum(w[1:num_surfaces]) : w[1:num_surfaces]
        F_raw = F_raw[1:num_surfaces,1:num_surfaces]
    else
        w = renorm ? w./minimum(w) : w
    end

    return DkAP(F_raw, w; k_dykstra=k_dykstra, max_iters=max_iters, verbose=verbose, nz_over_N=nz_over_N)
end

function AP_convergence_check(w::AbstractVector, num_surfaces::Int)
    num_volumes = num_surfaces == length(w) ? 0 : length(w) - num_surfaces
    if num_volumes == 0
        w_max = maximum(w)
        w_sum = sum(w)
        check_passed = w_max < 0.5 * w_sum
        if !check_passed
            error("Smoothing convergence check failed: max surface w ($w_max) ≥ half of total w ($(w_sum/2));\n
                   smoothing convergence not possible, refine the mesh (or trace more rays).")
        end
    end
end

function AP(F::AbstractMatrix, w::AbstractVector;
            max_iters::Int=1_000,
            verbose::Bool=true,
            M::Union{Nothing, AbstractMatrix}=nothing,
            nz_over_N::AbstractFloat=1.0,
            mult::AbstractFloat=1.0)

    # check if convergence is possible
    # count zero diagonal (zero self-view for planar surfaces)
    num_surfaces = sum([F[i,i] < 4*eps(Float64) for i in 1:length(w)])
    AP_convergence_check(w, num_surfaces)
    
    # Choose solver strategy
    if F isa SparseMatrixCSC
        use_sparse = true
    else
        use_sparse = false
        F = Matrix(F)
    end
    use_parallel = length(w) > 1000 && Threads.nthreads() > 1

    verbose && println("Matrix size: $(length(w))×$(length(w))")
    verbose && println("Strategy: $(use_parallel ? "Parallel" : "Serial"), $(use_sparse ? "sparse" : "dense")")

    if use_sparse
        F, tperm = symmetric_workspace(F)
        tmp    = similar(nonzeros(F))
        rowsum = zeros(eltype(F), length(w))
    elseif use_parallel
        steps_1_to_4! = all_steps_parallel!
    else
        steps_1_to_4! = all_steps!
    end

    # --- stopping machinery -------------------------------------------------
    N = length(w)
    target = 8 * eps(Float64)                          # best-case terminal condition
    guard  = 8 * sqrt(N / nz_over_N) * eps(Float64)    # floor-aware acceptance level
    switch = 4 * guard                                 # below this: check every iteration
    max_stride = 52                                    # cap on scheduled check stride: abs(log2(eps(Float64)))

    delta_prev = 100.0
    k_prev = 0
    k = 0
    k_next = 0
    c = 0
    flat = 0
    rho_est = 0.5
    converged = false
    floor_accepted = false

    delta_init = if use_sparse
        delta_perp(F, w; tperm=tperm, mode_indicator=:AP_sparse) # lower bound
    elseif M === nothing
        delta_perp(F, w; mode_indicator=:AP_sparse)              # lower bound
    else
        delta_perp(F, w; mode_indicator=:AP_dense)          # exact distance
    end
    delta_best = delta_init
    bounded = use_sparse || M === nothing
    if bounded
        verbose && println("Initial perpendicular distance to target manifold: $delta_init ≤ δ⟂ ≤ $(mult*delta_init)")
    else
        verbose && println("Initial perpendicular distance to target manifold: δ⟂ = $delta_init")
    end
    delta = delta_init

    while k < max_iters && delta > target
        if use_sparse
            all_steps_sparse!(F, w, tperm, tmp, rowsum)
        else
            steps_1_to_4!(F, w)
        end
        k += 1
        if k >= k_next
            delta = if use_sparse
                delta_perp(F, w; tperm=tperm, mode_indicator=:AP_sparse)
            else
                delta_perp(F, w, mode_indicator=:AP_dense)
            end
            c += 1

            # contraction-rate estimate over the elapsed stride
            if c >= 3
                rho_est = min(max((delta/delta_prev)^(1/max(k - k_prev, 1)), 0.5), 0.9999)
            end

            # flatness counter: three consecutive checks with <0.1% relative improvement
            if (c >= 3 && (delta_prev - delta) / max(delta_prev, eps(Float64)) < 1e-3) || (c >= 3 && delta >= delta_best * (1 - 1e-3))
                flat += 1
            else
                flat = 0
            end
            delta_best = min(delta_best, delta)

            # acceptance: tight stall-accept OR flatness (floor reached)
            if (rho_est > 0.99 && delta < guard) || flat >= 3
                floor_accepted = true
                verbose && println("Alternating projection (AP) contraction exhausted at iteration $k (floor δ_R ≈ $delta)")
                break
            end

            k_prev, delta_prev = k, delta

            # scheduling: rate-based stride above `switch`, per-iteration below;
            # stride capped so a clamped rho_est can never schedule past reach
            if delta > switch
                stride = max(1, ceil(Int, log(delta/target) / log(1/rho_est)))
                k_next = k + min(stride, max_stride)
            else
                k_next = k + 1
            end
            if bounded
                verbose && println("Iteration $k: perpendicular distance to target manifold: $delta ≤ δ⟂ ≤ $(mult*delta)")
            else
                verbose && println("Iteration $k: perpendicular distance to target manifold: δ⟂ = $delta")
            end
        end
    end
    converged = delta <= target || floor_accepted
    # ------------------------------------------------------------------------

    if converged
        if bounded
            verbose && println("Converged after $k iterations. Final perpendicular distance to manifold: $delta ≤ δ⟂ ≤ $(mult*delta)")
        else
            verbose && println("Converged after $k iterations. Final perpendicular distance to manifold: δ⟂ = $delta")
        end
    else
        if bounded
            @warn "Reached max_iters = $max_iters. Final perpendicular distance to manifold: $delta ≤ δ⟂ ≤ $(mult*delta)"
        else
            @warn "Reached max_iters = $max_iters. Final perpendicular distance to manifold: δ⟂ = $delta"
        end
    end

    if delta > delta_init && delta > sqrt(nz_over_N)
        @warn "Smoothing increased the perpendicular distance to the target manifold."
        @warn "Use unsmoothed mesh.F_raw instead of mesh.F_smooth."
    end

    return F
end