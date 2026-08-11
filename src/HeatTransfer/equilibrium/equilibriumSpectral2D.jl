function equilibriumSpectral2D!(rtm::RayTracingDomain2D, F_matrices::Union{AbstractMatrix, Vector{AbstractMatrix}}; 
                              max_iters::P=1000, convergence_tol=1e-12) where {P<:Integer}
    """
    Solve spectral radiative equilibrium using the optimal solver.
    
    Automatically chooses between:
    - Direct emission solver: for ε=1 everywhere AND σₛ=0 everywhere (orders of magnitude faster)
    - Full iterative solver: for general problems with scattering/reflection
    
    Set force_full=true to disable optimization and always use full solver.
    """
    G = Float64
    
    # Check if we can use the fast direct solver
    if rtm.spectral_mode == :spectral_uniform
        # :spectral_uniform routes to the direct solver, whose per-bin
        # reconstruction is only valid without reflection or scattering.
        # Auto-detection sets the mode correctly; this fires only if the
        # mode was set manually against the physics.
        b_check = get_b(rtm)
        if any(b_check .> 1e-12)
            error("""
            spectral_mode is :spectral_uniform, but the domain has reflecting
            surfaces or scattering volumes (max b = $(maximum(b_check))).
            The uniform-mode solver is not valid for this configuration and
            would produce incorrect per-bin results.

            The automatic detection selects the correct solver - avoid setting it manually.
            Remove the manual assignment, or set :spectral_variable.
            """)
        end
        println("=== Using DIRECT spectral solver ===")
        return equilibriumSpectral2D_direct!(rtm, F_matrices)
    elseif rtm.spectral_mode == :spectral_variable
        println("=== Using FULL spectral solver ===")
        return equilibriumSpectral2D_woodbury!(rtm, F_matrices;
                                          max_iters=max_iters,
                                          convergence_tol=convergence_tol)
    else
        error("spectral_mode must be :spectral_uniform or :spectral_variable but is $(rtm.spectral_mode)")
    end
end

function equilibriumSpectral2D_direct!(rtm::RayTracingDomain2D, F_matrices::AbstractMatrix) # single F
    """
    OPTIMIZED direct emission solver for non-scattering, non-reflecting problems.
    
    Key optimization: Instead of solving (n_bins+1)*N × n_bins*N block system for j,
    we solve N×N system for e directly, then recover j from emission fractions.
    
    This reduces computational complexity from O((n_bins*N)²) to O(N²).
    
    Requires: ε=1 everywhere (no reflection) AND σₛ=0 everywhere (no scattering)
    """
    G = Float64
    
    # Validate spectral setup at entry
    if isnothing(rtm.wavelength_band_limits)
        error("""
        Spectral solve requires wavelength band limits to be set.
        
        Please set mesh.wavelength_band_limits before calling steadyStateSpectral2D_direct!:
        
        Example, using logarithmic spacing:
            mesh.wavelength_band_limits = 10 .^ range(log10(0.0000001), log10(0.001), length=51)
        """)
    end
    
    if length(rtm.wavelength_band_limits) < 4
        error("wavelength_band_limits must have at least 4 values (defining 3 bins)")
    end
    
    if any(rtm.wavelength_band_limits .<= 0)
        error("wavelength_band_limits must all be positive (wavelengths > 0)")
    end

    if any(diff(rtm.wavelength_band_limits) .<= 0)
        error("wavelength_band_limits must be strictly increasing (no duplicates)")
    end

    # Get system dimensions
    surface_mapping, volume_mapping = rtm.surface_mapping, rtm.volume_mapping
    num_surfaces = length(surface_mapping)
    num_volumes = length(volume_mapping)
    total_elements = num_surfaces + num_volumes
    
    # Uniform spectral - single F matrix for all bins
    M = buildSystemMatrix(rtm, F_matrices; spectral_bin=1)
        
    # Set boundary conditions from mesh
    boundary, temperatures, emissive = setupBoundaryConditions(rtm, F_matrices)
    boundary = rtm.surfaces_only ? boundary[1:length(rtm.surface_mapping)] : boundary    
    b = get_b(rtm) # reflection-scattering

    # Iteration loop
    println("Starting spectral steady-state direct solve...")
    
    # the system simplifies
    j_tot = M \ boundary
    emissive = (I - Diagonal(b[:,1])*F_matrices')*j_tot
    
    # iterate to find the correct blackbody spectral fractions
    emitFrac = getBinsEmissionFractions(rtm, temperatures)
    temperatures = updateTemperaturesSpectral!(rtm, emissive, emitFrac)
    temperatures_prev = temperatures
    for iter = 1:10
        emitFrac = getBinsEmissionFractions(rtm, temperatures)
        temperatures = updateTemperaturesSpectral!(rtm, emissive, emitFrac)
        max_iteration_change = maximum(abs.(temperatures - temperatures_prev))
        if max_iteration_change < 1e-3
            break
        end
        temperatures_prev = temperatures
    end

    # Reconstruct final sol_j from final j_tot and emitFrac
    elements_tot = rtm.surfaces_only ? length(rtm.surface_mapping) : total_elements
    weightedFrac = getWeightedEmissionFractions(rtm, temperatures)
    sol_j = zeros(G, rtm.n_spectral_bins * elements_tot)
    for bin in 1:rtm.n_spectral_bins
        bin_start = (bin - 1) * elements_tot + 1
        bin_end = bin * elements_tot
        sol_j[bin_start:bin_end] = weightedFrac[1:elements_tot, bin] .* j_tot
    end
    
    # Write spectral results for each bin to mesh
    println("Writing spectral results to mesh...")
    for bin in 1:rtm.n_spectral_bins
        # Extract solution for this bin
        bin_start = (bin - 1) * elements_tot + 1
        bin_end = bin * elements_tot
        j_bin = sol_j[bin_start:bin_end]
        
        # Compute e, r, g_a from GERT matrices
        e_bin = (I - Diagonal(b[:,1])*F_matrices') * j_bin           # e = D * j
        r_bin = j_bin - e_bin                     # r = j - e = R' * j
        g_a_bin = j_bin - (I - F_matrices') * j_bin - r_bin  # g_a = A' * j = j - C*j - r # C
        
        # Write to mesh for this bin
        writeResultsToDomain!(rtm, j_bin, g_a_bin, r_bin; spectral_bin=bin)
    end

    # Update mesh with final results
    writeTemperaturesHeatSources!(rtm, temperatures)
    
    # Compute energy conservation error for each spectral bin
    rtm.energy_error = G.([sum((I - F_matrices') * sol_j[(i - 1) * elements_tot + 1:i * elements_tot]) 
                        for i in 1:rtm.n_spectral_bins])    
end

###### NEWTON-ACCELERATED WOODBURY SOLVER ######
# =============================================================================
# Sparse-capable Woodbury spectral solver for RayTraceHeatTransfer.jl
#
# Drop-in replacement for equilibriumSpectral2D_woodbury! plus support pieces.
# Dense path is byte-for-byte the existing algorithm; a sparse branch is added
# that mirrors the grey solver's philosophy: factor/iterate on the original
# sparsity pattern, never materialize squared or inverted operators.
#
# Key identity (D_k square, nonsingular since rho(R_k) < 1, Paper 1 Thm 4.7):
#     (D_k' D_k)^{-1} y  =  D_k^{-1} ( D_k^{-T} y )
# => one sparse LU of D_k per band replaces the dense Cholesky of D_k'D_k
#    (no pattern squaring, no condition-number squaring, no ridge needed).
# The Woodbury inner matrix  I + sum_k M_k (D_k'D_k)^{-1} M_k'  is SPD and is
# solved by CG with a matrix-free matvec, never materialized.
#
# ALSO REQUIRED (one-line fix elsewhere):
#   buildSystemMatrix(domain, F::Matrix{P})  -->  F::AbstractMatrix
#   (SparseMatrixCSC is not a Matrix; grey never hit this because it builds M
#    inline. For sparse F, `permutedims(F)` and `I - Diagonal(c)*permutedims(F)`
#    stay sparse, so the function body needs no other change.)
# =============================================================================

# ---- Matrix-free SPD operator shim for Krylov.cg ----------------------------
struct WoodburyInnerOp{Ff}
    N::Int
    f!::Ff          # f!(y, v) computes y = A*v in-place
end
Base.size(A::WoodburyInnerOp)         = (A.N, A.N)
Base.size(A::WoodburyInnerOp, ::Int)  = A.N
Base.eltype(::WoodburyInnerOp)        = Float64
LinearAlgebra.mul!(y::AbstractVector, A::WoodburyInnerOp, v::AbstractVector) = A.f!(y, v)

# ---- (I - J_red) restricted to rad-eq indices, matrix-free ------------------
struct NewtonKrylovOp{Ff}
    n_req::Int
    req::Vector{Int}
    jmv!::Ff                 # jred_matvec!(out_full, v_full)
    v_full::Vector{Float64}  # scratch, length N
    Jv_full::Vector{Float64} # scratch, length N
end
Base.size(A::NewtonKrylovOp)        = (A.n_req, A.n_req)
Base.size(A::NewtonKrylovOp, ::Int) = A.n_req
Base.eltype(::NewtonKrylovOp)       = Float64
function LinearAlgebra.mul!(y::AbstractVector, A::NewtonKrylovOp, x::AbstractVector)
    fill!(A.v_full, 0.0)
    A.v_full[A.req] .= x            # embed: zero on prescribed elements
    A.jmv!(A.Jv_full, A.v_full)     # full-space J·v
    y .= x .- A.Jv_full[A.req]      # extract: (I - J)|_req
    return y
end

"""
    equilibriumSpectral2D_woodbury!(rtm, F_matrices; max_iters, convergence_tol)

Woodbury spectral steady-state solver. Automatically dispatches on the storage
of `F_matrices`:

  * dense  `Matrix`          -> existing path (dense Cholesky of D'D, dense
                                inner Cholesky, threaded band loop), unchanged
  * `SparseMatrixCSC`        -> sparse path: per-band sparse LU of D_k,
                                matrix-free CG for the inner N x N system,
                                O(nnz)-memory throughout

Results agree between paths to the precision of the linear solves.
"""
function equilibriumSpectral2D_woodbury!(rtm::RayTracingDomain2D,
                                          F_matrices::Union{AbstractMatrix, Vector{<:AbstractMatrix}};
                                          max_iters::P=500,
                                          convergence_tol=0.001) where {P<:Integer}

    G = Float64

    # ---- Validate spectral setup (unchanged) ---------------------------------
    if isnothing(rtm.wavelength_band_limits)
        error("""
        Spectral solve requires wavelength band limits to be set.

        Please set mesh.wavelength_band_limits before calling
        equilibriumSpectral2D_woodbury!:

        Example, using logarithmic spacing:
            mesh.wavelength_band_limits = 10 .^ range(log10(0.0000001), log10(0.001), length=51)
        """)
    end
    if length(rtm.wavelength_band_limits) < 4
        error("wavelength_band_limits must have at least 4 values (defining 3 bins)")
    end
    if any(rtm.wavelength_band_limits .<= 0)
        error("wavelength_band_limits must all be positive (wavelengths > 0)")
    end
    if any(diff(rtm.wavelength_band_limits) .<= 0)
        error("wavelength_band_limits must be strictly increasing (no duplicates)")
    end

    # ---- System dimensions ---------------------------------------------------
    surface_mapping, volume_mapping = rtm.surface_mapping, rtm.volume_mapping
    num_surfaces = length(surface_mapping)
    num_volumes  = rtm.surfaces_only ? 0 : length(volume_mapping)
    N            = num_surfaces + num_volumes
    K            = rtm.n_spectral_bins

    # ---- Per-band system matrices --------------------------------------------
    M_matrices = Vector{AbstractMatrix}(undef, K)
    for bin in 1:K
        M_matrices[bin] = buildSystemMatrix(rtm, F_matrices[bin]; spectral_bin=bin)
    end

    b = get_b(rtm)
    use_sparse = F_matrices[1] isa SparseMatrixCSC

    # D_k kept explicitly in both paths (RHS assembly + result writing)
    D_matrices = Vector{AbstractMatrix}(undef, K)
    for k in 1:K
        Dk = I - Diagonal(b[:,k]) * F_matrices[k]'
        D_matrices[k] = use_sparse ? sparse(Dk) : Dk
    end

    # ---- Woodbury precomputation ---------------------------------------------
    println("==== Building and Factorizing Woodbury blocks ====")
    println("    Storage: $(use_sparse ? "sparse (LU per band + matrix-free CG inner)" :
                                          "dense (Cholesky per band + dense inner)")")

    local solve_inner            # Mu::Vector -> w::Vector
    local DtD_solve!             # (out, k, y, tmp) -> out = (Dk'Dk)^{-1} y
    local apply_DtDinvMt!        # (out, k, w, t1, t2) -> out = (Dk'Dk)^{-1} Mk' w
    # Dense-path objects (used by the dense predictor call at the end)
    local DtD_factors, DtDinv_Mt, inner_factor

    if use_sparse
        D_factors = Vector{Any}(undef, K)
        for k in 1:K
            D_factors[k] = lu(D_matrices[k]::SparseMatrixCSC)      # UMFPACK
        end

        DtD_solve! = (out, k, y, tmp) -> begin
            ldiv!(tmp, D_factors[k]', y)       # tmp = Dk^{-T} y
            ldiv!(out, D_factors[k],  tmp)     # out = Dk^{-1} tmp
            out
        end

        apply_DtDinvMt! = (out, k, w, t1, t2) -> begin
            mul!(t1, M_matrices[k]', w)
            DtD_solve!(out, k, t1, t2)
        end

        # inner matvec: y = v + sum_k Mk (Dk'Dk)^{-1} Mk' v   (serial over bands)
        cg_t1 = zeros(G, N); cg_t2 = zeros(G, N); cg_t3 = zeros(G, N)
        inner_mul! = (y, v) -> begin
            copyto!(y, v)
            for k in 1:K
                mul!(cg_t1, M_matrices[k]', v)
                DtD_solve!(cg_t2, k, cg_t1, cg_t3)
                mul!(y, M_matrices[k], cg_t2, one(G), one(G))
            end
            y
        end
        inner_op = WoodburyInnerOp(N, inner_mul!)
        solve_inner = Mu -> begin
            w, stats = Krylov.cg(inner_op, Mu; rtol = 1e-12, atol = 0.0)
            stats.solved || @warn "Woodbury inner CG did not fully converge" stats.niter
            w
        end
    else
        # ---- existing dense precomputation, unchanged ----
        DtD_factors = Vector{Cholesky{G, Matrix{G}}}(undef, K)
        DtDinv_Mt   = Vector{AbstractMatrix}(undef, K)
        ridge       = sqrt(eps(real(G))) * one(real(G))
        inner = Matrix{G}(I, N, N)
        for k in 1:K
            Dk   = D_matrices[k]
            DtDk = Symmetric(Dk' * Dk + ridge * I)
            chol = cholesky(DtDk)
            DtD_factors[k] = chol
            DtDinv_Mt[k]   = chol \ Matrix(M_matrices[k]')
            inner .+= M_matrices[k] * DtDinv_Mt[k]
        end
        inner_factor = cholesky(Symmetric(inner))

        DtD_solve!      = (out, k, y, tmp) -> ldiv!(out, DtD_factors[k], y)
        apply_DtDinvMt! = (out, k, w, t1, t2) -> mul!(out, DtDinv_Mt[k], w)
        solve_inner     = Mu -> inner_factor \ Mu
    end

    # ---- Boundary conditions and iteration state -----------------------------
    boundary, temperatures, emissive = setupBoundaryConditions(rtm, F_matrices)

    sol_j                   = zeros(G, K * N)
    previous_sol_j          = zeros(G, K * N)
    previous_previous_sol_j = zeros(G, K * N)
    # emitFrac                = getBinsEmissionFractions(rtm, temperatures)
    emitFrac                = getWeightedEmissionFractions(rtm, temperatures)

    u_list   = [zeros(G, N) for _ in 1:K]
    Mu_total = zeros(G, N)

    # Per-chunk buffers, hoisted OUT of the iteration loop.
    nthr    = Threads.nthreads() # use_sparse ? 1 :
    chunks  = collect(Iterators.partition(1:K, cld(K, nthr)))
    nchunks = length(chunks)
    Mu_partial = [zeros(G, N) for _ in 1:nchunks]
    e_k_buf    = [zeros(G, N) for _ in 1:nchunks]
    rhs_k_buf  = [zeros(G, N) for _ in 1:nchunks]
    tmp_k_buf  = [zeros(G, N) for _ in 1:nchunks]
    tmp2_k_buf = [zeros(G, N) for _ in 1:nchunks]

    # ---- Temperature fixed-point loop ---------------------------------------
    m_AA      = 5                       # Anderson depth
    ΔJ_hist   = Matrix{G}(undef, K*N, 0)
    Δf_hist   = Matrix{G}(undef, K*N, 0)
    f_prev    = zeros(G, K * N)
    j_in_prev = zeros(G, K * N)
    have_prev = false
    println("Starting spectral steady-state iteration (Woodbury)...")
    # One full sweep: given (emissive, emitFrac, temperatures) as INPUTS,
    # runs the Woodbury solve into sol_j, then updates emission state.
    # Returns (emissive_new, emitFrac_new, temperatures_new).
    function sweep!(emissive_in, emitFrac_in, temperatures_in, iter)
        h_top_l = rtm.surfaces_only ? boundary[1:length(rtm.surface_mapping)] : boundary

        for buf in Mu_partial; fill!(buf, zero(G)); end
        Threads.@sync for (t, chunk) in enumerate(chunks)
            Threads.@spawn begin
                e_k   = e_k_buf[t];  rhs_k = rhs_k_buf[t]
                tmp_k = tmp_k_buf[t]; Mu_t = Mu_partial[t]
                for k in chunk
                    Mk = M_matrices[k]; Dk = D_matrices[k]
                    @views e_k .= emissive_in .* emitFrac_in[:, k]
                    mul!(rhs_k, Mk', h_top_l)
                    mul!(rhs_k, Dk', e_k, one(G), one(G))
                    DtD_solve!(u_list[k], k, rhs_k, tmp_k)
                    mul!(Mu_t, Mk, u_list[k], one(G), one(G))
                end
            end
        end
        fill!(Mu_total, zero(G))
        for t in 1:nchunks; Mu_total .+= Mu_partial[t]; end
        w = solve_inner(Mu_total)
        Threads.@sync for (t, chunk) in enumerate(chunks)
            Threads.@spawn begin
                t1 = tmp_k_buf[t]; t2 = tmp2_k_buf[t]; corr = rhs_k_buf[t]
                for k in chunk
                    apply_DtDinvMt!(corr, k, w, t1, t2)
                    @views sol_j[(k-1)*N+1 : k*N] .= u_list[k] .- corr
                end
            end
        end

        em  = updateSpectralEmission!(rtm, iter, F_matrices, sol_j,
                                      emitFrac_in, temperatures_in, emissive_in)
        Tn  = updateTemperaturesSpectral!(rtm, em, getBinsEmissionFractions(rtm, temperatures_in))
        efn = getWeightedEmissionFractions(rtm, Tn)
        em  = rtm.surfaces_only ? em[1:length(rtm.surface_mapping)]     : em
        efn = rtm.surfaces_only ? efn[1:length(rtm.surface_mapping), :] : efn
        return em, efn, Tn, copy(sol_j)
    end
    newton_every  = 5                 # attempt Newton every 5th sweep
    J_newton      = nothing
    is_radeq_n    = falses(N)
    emissive_prev = copy(emissive)
    consec_rejects = 0; n_reassemble = 0
    for iter = 1:max_iters

        # one sweep step
        emissive, emitFrac, temperatures, sol_j = sweep!(emissive, emitFrac, temperatures, iter)

        # === Newton acceleration on total emissive =========================
        f_newton = emissive .- emissive_prev            # residual of the sweep
        if iter > 1 && iter % newton_every == 0
            if J_newton === nothing
                _, J_newton, is_radeq_n = predict_spectral_rate_dispatch(
                    rtm, F_matrices, M_matrices, D_matrices, b,
                    use_sparse,
                    use_sparse ? nothing : DtD_factors,
                    use_sparse ? nothing : DtDinv_Mt,
                    use_sparse ? nothing : inner_factor,
                    DtD_solve!, solve_inner,
                    sol_j, emissive, emitFrac, temperatures)
            end
            if J_newton !== nothing
                req = findall(is_radeq_n)
                δ = zeros(G, N)
                if J_newton isa AbstractMatrix
                    δ[req] .= (I - J_newton[req, req]) \ f_newton[req]
                elseif J_newton isa NewtonKrylovOp
                    δ_req, stats = Krylov.gmres(J_newton, f_newton[req];
                                                rtol = 1e-12, atol = 0.0)
                    stats.solved || @warn "Newton GMRES not fully converged" stats.niter
                    δ[req] .= δ_req
                end
                # probe: one sweep at e + λ·δ; returns residual norm and resulting state
                function ϕ_probe(λd)
                    e_t = max.(emissive .+ λd .* δ, zero(G))
                    T_t = updateTemperaturesSpectral!(rtm, e_t, getBinsEmissionFractions(rtm, temperatures))
                    ef_t = getWeightedEmissionFractions(rtm, T_t)
                    ef_t = rtm.surfaces_only ? ef_t[1:length(rtm.surface_mapping), :] : ef_t
                    em_v, ef_v, T_v, j_v = sweep!(copy(e_t), ef_t, T_t, iter)
                    return norm(em_v .- e_t), (copy(em_v), copy(ef_v), copy(T_v), j_v)
                end

                # --- capped line search: fixed candidates, plain step as floor ----
                ϕ0, s0 = ϕ_probe(0.0)                      # plain-step baseline
                ϕ_best, s_best, λ_best = ϕ0, s0, 0.0
                for λd in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0)
                    ϕt, st = ϕ_probe(λd)
                    if ϕt < ϕ_best
                        ϕ_best, s_best, λ_best = ϕt, st, λd
                    end
                end
                emissive, emitFrac, temperatures, j_b = s_best
                sol_j .= j_b
                if λ_best > 0
                    println("  Newton: accepted λ = $λ_best, ",
                            "residual $(round(ϕ_best/ϕ0, sigdigits=3))× plain step")
                    # after an acceptance: 
                    consec_rejects = 0
                else
                    println("  Newton: direction rejected (plain step kept)")
                    # after a rejection:  
                    consec_rejects += 1
                    if consec_rejects >= 3 # && n_reassemble < 3
                        J_newton = nothing; n_reassemble += 1; consec_rejects = 0
                    end
                end
            end
        end
        emissive_prev = copy(emissive)
        # ===================================================================

        f = sol_j .- previous_sol_j                 # residual of the raw solve
        convergence_error = norm(f) / norm(sol_j)
        previous_previous_sol_j .= previous_sol_j
        previous_sol_j .= sol_j

        if iter % 20 == 0
            println("Iteration $iter: convergence error = $convergence_error")
        end

        if iter > 1 && convergence_error < convergence_tol
            emissive, emitFrac, temperatures, sol_j = sweep!(emissive, emitFrac, temperatures, max_iters)
            emissive     = updateSpectralEmission!(rtm, iter, F_matrices, sol_j,
                                                   emitFrac, temperatures, emissive)
            temperatures = updateTemperaturesSpectral!(rtm, emissive, getBinsEmissionFractions(rtm, temperatures))
            println("Converged after $iter iterations, Final errors: convergence error = $convergence_error")
            break
        end

        if iter == max_iters
            println("Warning: Maximum iterations reached:")
            println("Final errors: convergence error = $convergence_error")
            emissive, emitFrac, temperatures, sol_j = sweep!(emissive, emitFrac, temperatures, max_iters)
            emissive     = updateSpectralEmission!(rtm, iter, F_matrices, sol_j,
                                                   emitFrac, temperatures, emissive)
            temperatures = updateTemperaturesSpectral!(rtm, emissive, getBinsEmissionFractions(rtm, temperatures))
        end
    end

    # ---- Write spectral results back to the domain ---------------------------
    # (adjoint spmv throughout: sparse-safe unchanged)
    println("Writing spectral results to mesh...")
    for bin in 1:K
        bin_start = (bin - 1) * N + 1
        bin_end   = bin * N
        j_bin     = sol_j[bin_start:bin_end]

        e_bin   = D_matrices[bin] * j_bin
        r_bin   = j_bin .- e_bin
        g_a_bin = j_bin .- (I - F_matrices[bin]') * j_bin .- r_bin

        writeResultsToDomain!(rtm, j_bin, g_a_bin, r_bin; spectral_bin=bin)
    end

    writeTemperaturesHeatSources!(rtm, temperatures)

    rtm.energy_error = G.([sum((I - F_matrices[i]') *
                               sol_j[(i - 1) * N + 1 : i * N]) for i in 1:K])

    return rtm
end


"""
    predict_spectral_convergence_rate(rtm, D_matrices, M_matrices,
                                       DtD_factors, DtDinv_Mt, inner_factor,
                                       sol_j, emissive, emitFrac, temperatures;
                                       use_reduced=true)

Compute the predicted asymptotic spectral radius ρ(J) of the spectral Woodbury
fixed-point iteration at the current iterate.

The fixed-point Jacobian is

    J  =  H^{-1} D̃ᵀ · 𝓔_J · D_sum                            (K·N × K·N)

with reduced form sharing the same nonzero spectrum

    J_red  =  D_sum · H^{-1} D̃ᵀ · 𝓔_J                         (N × N)

where:
  • D_sum  = [D_1 | D_2 | … | D_K]          (N × K·N)
  • D̃     = blkdiag(D_k)                   (K·N × K·N)
  • H      = D̃ᵀ D̃ + M̃ᵀ M̃                  (already factored in Woodbury form)
  • 𝓔_J    = stack of K blocks; block k = F_k + diag(e) · Φ_k · J_T   (N × N)
      F_k = diag(emitFrac[:,k])
      Φ_k = diag(∂f_k/∂T at current T), zero on prescribed-T elements
      J_T = diag(T_i / (4 e_i)), zero on prescribed-T elements
      e   = current `emissive` vector  (per-element total emissive power)

The action  H^{-1} D̃ᵀ v   is computed via the SAME Woodbury blocks used by the
solver — DtD_factors, DtDinv_Mt, inner_factor — at one Woodbury solve per
matvec. The dominant eigenvalue is found by power iteration on J_red, which is
N × N and dense and small enough that direct eigensolve is also viable.

Pass `use_reduced=false` to assemble and eigensolve the full K·N × K·N J
instead (more expensive; mainly for cross-validation).

Returns: `rho_pred::Float64` — predicted asymptotic ρ.
"""
function predict_spectral_convergence_rate(rtm::RayTracingDomain2D,
        F_matrices::Vector{AbstractMatrix},
        M_matrices::Vector{AbstractMatrix},
        DtD_factors::Vector{<:Cholesky},
        DtDinv_Mt::Vector{AbstractMatrix},
        inner_factor::Cholesky,
        sol_j::Vector{G},
        emissive::Vector{G},
        emitFrac::AbstractMatrix,
        temperatures::Vector{G};
        use_reduced::Bool=true) where {G}

    N = length(emissive)
    K = rtm.n_spectral_bins
    b = get_b(rtm)

    # --- 1. Determine which elements are radiative-equilibrium (Jacobian-active) ---
    #     prescribed-T elements have ∂e/∂j = 0 ⇒ zero row in J_T
    is_radeq = falses(N)
    num_surfaces = length(rtm.surface_mapping)
    for ((ci, fi, wi), si) in rtm.surface_mapping
        face = rtm.fine_mesh[ci][fi]
        is_radeq[si] = face.T_in_w[wi] < -0.1
    end
    if !rtm.surfaces_only
        for ((ci, fi), vi) in rtm.volume_mapping
            face = rtm.fine_mesh[ci][fi]
            is_radeq[num_surfaces + vi] = face.T_in_g < -0.1
        end
    end

    # --- 2. J_T diagonal: T_i / (4 e_i) on rad-eq elements, 0 elsewhere ---
    J_T_diag = G.(ifelse.(is_radeq, temperatures ./ (4 .* emissive), zero(G)))

    # --- 3. Φ_k = diag(∂f_k/∂T) at current T (one column per band) ---
    #     ∂f_k/∂T = λ_k · F'(λ_k T) − λ_{k−1} · F'(λ_{k−1} T)
    #     with F'(λT) = dF_blackbody/d(λT) evaluated at λT.
    Phi = zeros(G, N, K)
    for i in 1:N
        T_i = temperatures[i]
        T_i <= 0 && continue
        lam = rtm.wavelength_band_limits
        for k in 1:K
            dF_lo = (k == 1) ? zero(G) : lam[k]   * dF_blackbody_dlambdaT(lam[k]   * T_i)
            dF_hi = (k == K) ? zero(G) : lam[k+1] * dF_blackbody_dlambdaT(lam[k+1] * T_i)
            Phi[i, k] = dF_hi - dF_lo
        end
    end

    # --- 4. Define the matvec  v ∈ ℝᴺ ↦ J_red · v ∈ ℝᴺ ---
    #     a) build 𝓔_J · v band-by-band:
    #          (𝓔_J v)_k = (emitFrac[:,k] .+ emissive .* Phi[:,k] .* J_T_diag) .* v
    #     b) RHS_k = D_k' (𝓔_J v)_k
    #     c) Woodbury solve to get j_vec (K·N)
    #     d) output = D_sum · j_vec = Σ_k D_k · j_vec_k
    function jred_matvec(v::AbstractVector)
        Mu_total = zeros(G, N)
        u_list   = [zeros(G, N) for _ in 1:K]
        for k in 1:K
            ek_v       = (emitFrac[:, k] .+ emissive .* Phi[:, k] .* J_T_diag) .* v
            rhs_k      = (I - Diagonal(b[:,k])*F_matrices[k]') * ek_v
            u_list[k] .= DtD_factors[k] \ rhs_k
            Mu_total .+= M_matrices[k] * u_list[k]
        end
        w = inner_factor \ Mu_total
        out = zeros(G, N)
        for k in 1:K
            jk   = u_list[k] .- DtDinv_Mt[k] * w
            out .+= (I - Diagonal(b[:,k])*F_matrices[k]') * jk
        end
        return out
    end

    # --- 5. Get ρ ---
    if use_reduced
        # Materialize J_red (N × N, small) by applying matvec to N basis vectors,
        # then dense eigensolve. For larger N, swap in KrylovKit.eigsolve.
        J_red = zeros(G, N, N)
        e_basis = zeros(G, N)
        for j in 1:N
            fill!(e_basis, zero(G)); e_basis[j] = one(G)
            J_red[:, j] .= jred_matvec(e_basis)
        end
        # return Float64(maximum(abs.(eigvals(J_red))))
        return Float64(-1), J_red, is_radeq # Float64(maximum(abs.(eigvals(J_red))))
    else
        # Full K·N × K·N assembly (for cross-validation only; do not use for big problems)
        Dtil = zeros(G, K*N, K*N)
        for k in 1:K
            Dtil[(k-1)*N+1:k*N, (k-1)*N+1:k*N] .= I - Diagonal(b[:,k])*F_matrices[k]'
        end
        Dsum = hcat([I - Diagonal(b[:,k])*F_matrices[k]' for k in 1:K]...)
        EJ   = zeros(G, K*N, N)
        for k in 1:K
            EJ[(k-1)*N+1:k*N, :] .= Diagonal(emitFrac[:, k]) .+
                                     Diagonal(emissive) * Diagonal(Phi[:, k]) * Diagonal(J_T_diag)
        end
        Mtil = hcat(M_matrices...)
        H    = Dtil' * Dtil + Mtil' * Mtil
        J    = H \ (Dtil' * EJ * Dsum)
        return Float64(-1) # Float64(maximum(abs.(eigvals(J))))
    end
end

# =============================================================================
# Rate predictor dispatch: dense path calls the existing implementation
# unchanged; sparse/large path runs power iteration on the SAME matrix-free
# jred_matvec (never materializes J_red).
# =============================================================================
function predict_spectral_rate_dispatch(rtm, F_matrices, M_matrices, D_matrices, b,
                                        use_sparse,
                                        DtD_factors, DtDinv_Mt, inner_factor,
                                        DtD_solve!, solve_inner,
                                        sol_j, emissive, emitFrac, temperatures)
    G = Float64
    N = length(emissive)
    K = rtm.n_spectral_bins

    if !use_sparse && N <= 2000
        return predict_spectral_convergence_rate(
            rtm, F_matrices, M_matrices,
            DtD_factors, DtDinv_Mt, inner_factor,
            sol_j, emissive, emitFrac, temperatures)
    end

    # ---- shared setup (mirrors predict_spectral_convergence_rate steps 1-3) --
    is_radeq = falses(N)
    num_surfaces = length(rtm.surface_mapping)
    for ((ci, fi, wi), si) in rtm.surface_mapping
        face = rtm.fine_mesh[ci][fi]
        is_radeq[si] = face.T_in_w[wi] < -0.1
    end
    if !rtm.surfaces_only
        for ((ci, fi), vi) in rtm.volume_mapping
            face = rtm.fine_mesh[ci][fi]
            is_radeq[num_surfaces + vi] = face.T_in_g < -0.1
        end
    end
    J_T_diag = G.(ifelse.(is_radeq, temperatures ./ (4 .* emissive), zero(G)))

    Phi = zeros(G, N, K)
    for i in 1:N
        T_i = temperatures[i]
        T_i <= 0 && continue
        lam = rtm.wavelength_band_limits
        for k in 1:K
            dF_lo = (k == 1) ? zero(G) : lam[k]   * dF_blackbody_dlambdaT(lam[k]   * T_i)
            dF_hi = (k == K) ? zero(G) : lam[k+1] * dF_blackbody_dlambdaT(lam[k+1] * T_i)
            Phi[i, k] = dF_hi - dF_lo
        end
    end

    # ---- matrix-free J_red matvec via the solver's own hooks -----------------
    t1 = zeros(G, N); t2 = zeros(G, N); t3 = zeros(G, N)
    u_list = [zeros(G, N) for _ in 1:K]
    function jred_matvec!(out, v)
        Mu_total = zeros(G, N)
        for k in 1:K
            ek_v = (view(emitFrac, :, k) .+ emissive .* view(Phi, :, k) .* J_T_diag) .* v
            mul!(t1, D_matrices[k]', ek_v)          # rhs_k = Dk' (E_J v)_k
            DtD_solve!(u_list[k], k, t1, t2)
            mul!(Mu_total, M_matrices[k], u_list[k], one(G), one(G))
        end
        w = solve_inner(Mu_total)
        fill!(out, zero(G))
        for k in 1:K
            mul!(t1, M_matrices[k]', w)
            DtD_solve!(t2, k, t1, t3)
            t2 .= u_list[k] .- t2                   # j_k
            mul!(out, D_matrices[k], t2, one(G), one(G))   # out += Dk j_k
        end
        out
    end

    # ---- materialize J_red column-by-column via the matrix-free matvec ------
    J_red = zeros(G, N, N)
    e_basis = zeros(G, N)
    col = zeros(G, N)
    for j in 1:N
        fill!(e_basis, zero(G)); e_basis[j] = one(G)
        jred_matvec!(col, e_basis)
        J_red[:, j] .= col
    end
    rho = -1.0 # maximum(abs.(eigvals(J_red)))
    # ---- wrap the matvec as a restricted operator for Newton–Krylov ---------
    req_idx = findall(is_radeq)
    jred_op = NewtonKrylovOp(length(req_idx), req_idx, jred_matvec!,
                             zeros(G, N), zeros(G, N))
    return Float64(rho), jred_op, is_radeq

end

predict_rho(args...; kwargs...) = begin
    J, _ = predict_spectral_convergence_rate(args...; kwargs...)
    Float64(maximum(abs.(eigvals(J))))
end

"""
    dF_blackbody_dlambdaT(x)

Derivative dF_(0→λT)/d(λT) of the fractional blackbody function evaluated at
argument x = λT. Equals the normalized Planck spectrum, i.e. the integrand of
the F(λT) integral. Uses the standard exponentially-convergent series.
"""
function dF_blackbody_dlambdaT(λT::Real)
    λT <= 0 && return zero(λT)
    C2 = 1.4387769e-2        # m·K — must match what emitFracBlackBodySpectrum uses
    x  = C2 / λT
    x > 200 && return zero(λT)
    # d/d(λT) of (15/π⁴) Σ exp(-nx)/n · (x³ + 3x²/n + 6x/n² + 6/n³)
    # Use the cleaner closed form: dF/d(λT) = (15/π⁴) · C2⁴ / (λT)⁵ · 1/(exp(x) − 1)
    # (this is the Planck spectrum integrand in the F-of-λT variable).
    denom = expm1(x)
    return (15 / π^4) * C2^4 / (λT^5 * denom)
end

### ORIGINAL FIXED-POINT SOLVER (kept for reference) ###
# function equilibriumSpectral2D_full!(rtm::RayTracingDomain2D, F_matrices::Union{AbstractMatrix, Vector{AbstractMatrix}}; 
#                                     max_iters::P=500,
#                                     convergence_tol=0.001) where {P<:Integer}
#     """
#     Solve spectral radiative equilibrium using the iterative method from your example.
#     F_matrices can be either a single matrix (grey/uniform) or vector of matrices (variable spectral)
#     """
#     G = Float64
    
#     # Validate spectral setup at entry
#     if isnothing(rtm.wavelength_band_limits)
#         error("""
#         Spectral solve requires wavelength band limits to be set.
        
#         Please set mesh.wavelength_band_limits before calling steadyStateSpectral2D!:
        
#         Example, using logarithmic spacing:
#             mesh.wavelength_band_limits = 10 .^ range(log10(0.0000001), log10(0.001), length=51)
#         """)
#     end
    
#     # Additional validation
#     if length(rtm.wavelength_band_limits) < 4
#         error("wavelength_band_limits must have at least 4 values (defining 3 bins)")
#     end
    
#     if any(rtm.wavelength_band_limits .<= 0)
#         error("wavelength_band_limits must all be positive (wavelengths > 0)")
#     end

#     if any(diff(rtm.wavelength_band_limits) .<= 0)
#         error("wavelength_band_limits must be strictly increasing (no duplicates)")
#     end

#     # Get system matrices
#     surface_mapping, volume_mapping = rtm.surface_mapping, rtm.volume_mapping
#     num_surfaces = length(surface_mapping)
#     num_volumes = length(volume_mapping)
#     total_elements = num_surfaces + num_volumes
    
#     # Build system matrix based on spectral mode
#     b = get_b(rtm) # get reflection-scattering coefficients
#     M_matrices = Vector{AbstractMatrix}()
#     # Variable spectral - use different F matrix for each bin
#     for bin in 1:rtm.n_spectral_bins
#         M = buildSystemMatrix(rtm, F_matrices[bin]; spectral_bin=bin)
#         push!(M_matrices, M)
#     end
    
#     # Build block matrix structure
#     println("==== Building and Factorizing Block matrix ====")
#     if F_matrices isa SparseMatrixCSC || (F_matrices isa Vector && F_matrices[1] isa SparseMatrixCSC)
#         block_matrix = spzeros((rtm.n_spectral_bins + 1) * total_elements, rtm.n_spectral_bins * total_elements)
#     else
#         block_matrix = zeros(G, (rtm.n_spectral_bins + 1) * total_elements, rtm.n_spectral_bins * total_elements)
#     end

#     for i = 1:(rtm.n_spectral_bins + 1)  # block indices
#         for j = 1:rtm.n_spectral_bins      # block indices
#             row_start = (i - 1) * total_elements + 1
#             row_end = i * total_elements
#             col_start = (j - 1) * total_elements + 1
#             col_end = j * total_elements
            
#             if i == 1
#                 block_matrix[row_start:row_end, col_start:col_end] = M_matrices[j]
#             elseif i == j + 1
#                 block_matrix[row_start:row_end, col_start:col_end] = I - Diagonal(b[:,j])*F_matrices[j]'
#             end
#             # else: leave as sparse zeros
#         end
#     end
#     # Outside iteration loop - factorize once
#     Factorization = qr(block_matrix)
    
#     # Set boundary conditions from mesh
#     boundary, temperatures, emissive = setupBoundaryConditions(rtm, F_matrices)

#     # Iteration loop
#     sol_j = zeros(G, rtm.n_spectral_bins * total_elements)
#     previous_sol_j = zeros(G, rtm.n_spectral_bins * total_elements)
#     emitFrac = getBinsEmissionFractions(rtm, temperatures)

#     println("Starting spectral steady-state iteration...")
#     for iter = 1:max_iters

#         # calculate emissive powers and temperatures
#         emissive = updateSpectralEmission!(rtm, iter, F_matrices, sol_j, emitFrac, temperatures, emissive)
#         temperatures = updateTemperaturesSpectral!(rtm, emissive, emitFrac)
#         emitFrac = getBinsEmissionFractions(rtm, temperatures)
        
#         # build and solve least-squares problem
#         emissive_pow_vec = reduce(vcat, [emissive for _ in 1:rtm.n_spectral_bins])
#         b_e_matrix = Diagonal(reduce(vcat, [boundary; emissive_pow_vec]))
#         sol_j .= Factorization \ (b_e_matrix*[ones(G, total_elements); emitFrac[:]])

#         # Check convergence
#         convergence_error = maximum(abs.(sol_j - previous_sol_j)) / maximum(abs.(sol_j))
#         previous_sol_j .= sol_j
#         if iter % 20 == 0
#             println("Iteration $iter: convergence error = $convergence_error")
#         end
        
#         if iter > 1 && convergence_error < convergence_tol
#             println("Converged after $iter iterations")
#             emissive = updateSpectralEmission!(rtm, iter, F_matrices, sol_j, emitFrac, temperatures, emissive)
#             temperatures = updateTemperaturesSpectral!(rtm, emissive, emitFrac)
#             break
#         end
        
#         if iter == max_iters
#             println("Warning: Maximum iterations reached:")
#             println("Final errors: convergence error = $convergence_error")
#             emissive = updateSpectralEmission!(rtm, iter, F_matrices, sol_j, emitFrac, temperatures, emissive)
#             temperatures = updateTemperaturesSpectral!(rtm, emissive, emitFrac)
#         end
#     end
    
#     # Write spectral results for each bin to mesh
#     println("Writing spectral results to mesh...")
#     for bin in 1:rtm.n_spectral_bins
#         # Extract solution for this bin
#         bin_start = (bin - 1) * total_elements + 1
#         bin_end = bin * total_elements
#         j_bin = sol_j[bin_start:bin_end]
        
#         # Compute e, r, g_a from GERT matrices
#         e_bin = (I - Diagonal(b)*F_matrices[bin]') * j_bin           # e = D * j
#         r_bin = j_bin - e_bin                     # r = j - e = R' * j
#         g_a_bin = j_bin - (I - F_matrices[bin]') * j_bin - r_bin  # g_a = A' * j = j - C*j - r # C_matrices[bin]
        
#         # Write to mesh for this bin
#         writeResultsToDomain!(rtm, j_bin, g_a_bin, r_bin; spectral_bin=bin)
#     end

#     # Update mesh with final results
#     writeTemperaturesHeatSources!(rtm, temperatures)

#     # Compute energy conservation error for each spectral bin
#     rtm.energy_error = G.([sum((I - F_matrices[i]') * sol_j[(i - 1) * total_elements + 1:i * total_elements]) 
#                             for i in 1:rtm.n_spectral_bins]) # C_matrices[i]
    
# end

###### ORIGINAL WOODBURY SOLVER (kept for reference) ######
# """
#     equilibriumSpectral2D_woodbury!(rtm, F_matrices; max_iters, convergence_tol)

# Drop-in alternative to `equilibriumSpectral2D_full!` that replaces the
# (K+1)*N × K*N QR factorization with a Woodbury decomposition.

# The block matrix has the structure

#     BM = [ M_1  M_2  ...  M_K ]    <- N rows (top stripe, couples all bands)
#          [ D_1                ]
#          [      D_2           ]    <- K*N rows (block-diagonal per band)
#          [          ...       ]
#          [               D_K  ]

# Decomposing BM = [M; D_blk] with M = [M_1|...|M_K] (size N x K*N) and
# D_blk = block_diag(D_k) (size K*N x K*N), the least-squares normal equations are

#     (D_blk^T D_blk + M^T M) j = M^T h_top + D_blk^T e_blk .

# Because M has only N rows, M^T M is rank-N, and Woodbury gives

#     (D_blk^T D_blk + M^T M)^{-1} = (D_blk^T D_blk)^{-1}
#         - (D_blk^T D_blk)^{-1} M^T * [I_N + M (D_blk^T D_blk)^{-1} M^T]^{-1} * M (D_blk^T D_blk)^{-1}.

# The block-diagonal (D_k^T D_k)^{-1} factors and the inner N x N inverse are
# pre-computed once before the temperature-iteration loop. The per-iteration cost
# drops from O(K^3 N^3) (full QR) to O(K N^2), i.e. a K^2 speedup.

# This solver produces results that agree with `equilibriumSpectral2D_full!` to
# the precision of the linear solve (~1e-10 relative for typical problems), with
# the same machine-precision per-band energy conservation when b is uniform
# across elements within each band.

# # Arguments
# - `rtm::RayTracingDomain2D`         : domain to solve on
# - `F_matrices`                      : vector of K F-matrices (one per band)
# - `max_iters::Integer=500`     : cap on fixed-point iterations
# - `convergence_tol=0.001`           : same criterion as the full solver
#                                       max|sol_j - prev_sol_j| / max|sol_j| < tol
# """
# function equilibriumSpectral2D_woodbury!(rtm::RayTracingDomain2D,
#                                           F_matrices::Union{AbstractMatrix, Vector{AbstractMatrix}};
#                                           max_iters::P=500,
#                                           convergence_tol=0.001) where {P<:Integer}

#     G = Float64

#     # ---- Validate spectral setup (mirrors equilibriumSpectral2D_full!) -------
#     if isnothing(rtm.wavelength_band_limits)
#         error("""
#         Spectral solve requires wavelength band limits to be set.

#         Please set mesh.wavelength_band_limits before calling
#         equilibriumSpectral2D_woodbury!:

#         Example, using logarithmic spacing:
#             mesh.wavelength_band_limits = 10 .^ range(log10(0.0000001), log10(0.001), length=51)
#         """)
#     end
#     if length(rtm.wavelength_band_limits) < 4
#         error("wavelength_band_limits must have at least 4 values (defining 3 bins)")
#     end
#     if any(rtm.wavelength_band_limits .<= 0)
#         error("wavelength_band_limits must all be positive (wavelengths > 0)")
#     end
#     if any(diff(rtm.wavelength_band_limits) .<= 0)
#         error("wavelength_band_limits must be strictly increasing (no duplicates)")
#     end

#     # ---- System dimensions ---------------------------------------------------
#     surface_mapping, volume_mapping = rtm.surface_mapping, rtm.volume_mapping
#     num_surfaces = length(surface_mapping)
#     num_volumes  = rtm.surfaces_only ? 0 : length(volume_mapping)
#     N            = num_surfaces + num_volumes      # total elements per band
#     K            = rtm.n_spectral_bins             # number of bands

#     # ---- Build per-band system matrices --------------------------------------
#     # D_matrices = Vector{AbstractMatrix}(undef, K)
#     # C_matrices = Vector{Matrix{G}}(undef, K)
#     M_matrices = Vector{AbstractMatrix}(undef, K)
#     for bin in 1:K
#         M = buildSystemMatrix(rtm, F_matrices[bin]; spectral_bin=bin) # C, D
#         # C_matrices[bin] = C
#         # D_matrices[bin] = D
#         M_matrices[bin] = M
#     end

#     # ---- Woodbury precomputation ---------------------------------------------
#     println("==== Building and Factorizing Woodbury blocks ====")

#     # For each band k:
#     #   - Cholesky-factor (D_k^T D_k)
#     #   - Pre-compute (D_k^T D_k)^{-1} M_k^T   (size N x N)
#     # These are reused at every iteration as matvecs.
#     #
#     # A small ridge (sqrt(eps) * I) is added to (D_k^T D_k) to handle
#     # near-singular cases robustly. In the well-conditioned regime
#     # (rho(R_k) bounded away from 1, per Paper 1 Theorem 4.7) this is
#     # numerically inactive.

#     DtD_factors = Vector{Cholesky{G, Matrix{G}}}(undef, K)
#     DtDinv_Mt   = Vector{AbstractMatrix}(undef, K)
#     ridge       = sqrt(eps(real(G))) * one(real(G))

#     inner = Matrix{G}(I, N, N)  # will accumulate I + sum_k M_k (D_k^T D_k)^{-1} M_k^T
#     b = get_b(rtm) # reflection-scattering

#     for k in 1:K
#         Dk     = I - Diagonal(b[:,k])*F_matrices[k]'
#         Mk     = M_matrices[k]
#         DtDk   = Symmetric(Dk' * Dk + ridge * I)
#         chol   = cholesky(DtDk)
#         DtD_factors[k] = chol
#         # (D_k^T D_k)^{-1} M_k^T  — N x N
#         DtDinv_Mt[k]   = chol \ Matrix(Mk')
#         # accumulate inner matrix
#         inner .+= Mk * DtDinv_Mt[k]
#     end

#     # Inner inverse: symmetric positive definite (I + sum of PSD terms)
#     inner_factor = cholesky(Symmetric(inner))

#     # ---- Boundary conditions and iteration state -----------------------------
#     boundary, temperatures, emissive = setupBoundaryConditions(rtm, F_matrices)

#     sol_j          = zeros(G, K * N)
#     previous_sol_j = zeros(G, K * N)
#     previous_previous_sol_j = zeros(G, K * N)
#     emitFrac       = getBinsEmissionFractions(rtm, temperatures)

#     # Per-band workspace buffers (avoid allocation in the loop)
#     u_list   = [zeros(G, N) for _ in 1:K]   # u_k = (D_k^T D_k)^{-1} (M_k^T h_top + D_k^T e_k)
#     Mu_total = zeros(G, N)                  # sum_k M_k u_k

#     # ---- Temperature fixed-point loop ---------------------------------------
#     println("Starting spectral steady-state iteration (Woodbury)...")
#     convergence_rate = Float64[]
#     for iter = 1:max_iters

#         # Update emissive power and temperature from previous j (same hooks as
#         # equilibriumSpectral2D_full!). On the first iteration sol_j is zero,
#         # which matches the full solver's behaviour.
#         emissive     = updateSpectralEmission!(rtm, iter, F_matrices, sol_j,
#                                                emitFrac, temperatures, emissive)
#         temperatures = updateTemperaturesSpectral!(rtm, emissive, emitFrac)
#         emitFrac     = getBinsEmissionFractions(rtm, temperatures)
#         emissive = rtm.surfaces_only ? emissive[1:length(rtm.surface_mapping)] : emissive
#         emitFrac = rtm.surfaces_only ? emitFrac[1:length(rtm.surface_mapping),:] : emitFrac

#         # Right-hand side of the least-squares system, expressed via the same
#         # scaling vector as the full solver:
#         #   h = b_e_matrix * [ones(N); emitFrac[:]]
#         # with b_e_matrix = Diagonal([boundary; repeat(emissive, K)])
#         # The top N entries of h are h_top = boundary .* 1 = boundary.
#         # Band-k entries are h_k = emissive .* emitFrac[:, k] = e_k.
#         h_top = rtm.surfaces_only ? boundary[1:length(rtm.surface_mapping)] : boundary
#         # e_k = emissive .* emitFrac[:, k] ; constructed inline below

#         # Hoisted ONCE outside the iteration loop (alongside u_list, Mu_total):
#         nthr = Threads.nthreads()
#         chunks = collect(Iterators.partition(1:K, cld(K, nthr)))
#         nchunks = length(chunks)   # may be < nthr if K < nthr; that's fine
#         Mu_partial = [zeros(G, N) for _ in 1:nchunks]
#         e_k_buf    = [zeros(G, N) for _ in 1:nchunks]
#         rhs_k_buf  = [zeros(G, N) for _ in 1:nchunks]

#         # === Woodbury solve =================================================
#         # Step 1 (per band):  u_k = (D_k^T D_k)^{-1} (M_k^T h_top + D_k^T e_k)
#         for buf in Mu_partial; fill!(buf, zero(G)); end

#         Threads.@sync for (t, chunk) in enumerate(chunks)
#             Threads.@spawn begin
#                 e_k   = e_k_buf[t]
#                 rhs_k = rhs_k_buf[t]
#                 Mu_t  = Mu_partial[t]
#                 for k in chunk
#                     Mk = M_matrices[k]
#                     Dk =  I - Diagonal(b[:,k])*F_matrices[k]'
#                     @views e_k .= emissive .* emitFrac[:, k]
#                     mul!(rhs_k, Mk', h_top)                # rhs_k = Mk' h_top
#                     mul!(rhs_k, Dk', e_k, one(G), one(G))  # rhs_k += Dk' e_k
#                     ldiv!(u_list[k], DtD_factors[k], rhs_k)
#                     mul!(Mu_t, Mk, u_list[k], one(G), one(G))  # Mu_t += Mk * u_list[k]
#                 end
#             end
#         end

#         # Reduce per-chunk partials into Mu_total (cheap: nchunks adds of length N)
#         fill!(Mu_total, zero(G))
#         for t in 1:nchunks
#             Mu_total .+= Mu_partial[t]
#         end

#         # Step 2: w = inner^{-1} * (sum_k M_k u_k)
#         w = inner_factor \ Mu_total

#         # Step 3 (per band):  j_k = u_k - (D_k^T D_k)^{-1} M_k^T * w
#         Threads.@sync for chunk in chunks
#             Threads.@spawn begin
#                 for k in chunk
#                     bin_start = (k - 1) * N + 1
#                     bin_end   = k * N
#                     @views sol_j[bin_start:bin_end] .= u_list[k] .- DtDinv_Mt[k] * w
#                 end
#             end
#         end
#         # ===================================================================

#         # # === Woodbury solve =================================================
#         # # Step 1 (per band):  u_k = (D_k^T D_k)^{-1} (M_k^T h_top + D_k^T e_k)
#         # fill!(Mu_total, zero(G))
#         # for k in 1:K
#         #     Mk  = M_matrices[k]
#         #     Dk  = D_matrices[k]
#         #     e_k = emissive .* @view emitFrac[:, k]
#         #     rhs_k        = Mk' * h_top + Dk' * e_k
#         #     u_list[k]   .= DtD_factors[k] \ rhs_k
#         #     Mu_total   .+= Mk * u_list[k]
#         # end

#         # # Step 2: w = inner^{-1} * (sum_k M_k u_k)
#         # w = inner_factor \ Mu_total

#         # # Step 3 (per band):  j_k = u_k - (D_k^T D_k)^{-1} M_k^T * w
#         # for k in 1:K
#         #     bin_start = (k - 1) * N + 1
#         #     bin_end   = k * N
#         #     sol_j[bin_start:bin_end] .= u_list[k] .- DtDinv_Mt[k] * w
#         # end
#         # # ===================================================================

#         # Convergence check — identical criterion to the full solver
#         if iter > 2
#             push!(convergence_rate, norm(sol_j .- previous_sol_j)/norm(previous_sol_j .-previous_previous_sol_j))
#         end
#         convergence_error = maximum(abs.(sol_j - previous_sol_j)) / maximum(abs.(sol_j))
#         previous_previous_sol_j .= previous_sol_j
#         previous_sol_j .= sol_j

#         if iter % 20 == 0
#             println("Iteration $iter: convergence error = $convergence_error")
#         end

#         if iter > 1 && convergence_error < convergence_tol
#             emissive     = updateSpectralEmission!(rtm, iter, F_matrices, sol_j,
#                                                    emitFrac, temperatures, emissive)
#             temperatures = updateTemperaturesSpectral!(rtm, emissive, emitFrac)
#             mean_rate = sum(convergence_rate)/length(convergence_rate)
#             rho_pred = predict_spectral_convergence_rate(
#                 rtm, F_matrices, M_matrices,
#                 DtD_factors, DtDinv_Mt, inner_factor,
#                 sol_j, emissive, emitFrac, temperatures)
#             println("Converged after $iter iterations, at rate $(round(mean_rate; sigdigits=10))\n
#                     Predicted rate $(round(rho_pred; sigdigits=10)), ratio: $(round(mean_rate/rho_pred; sigdigits=10))")
#             break
#         end

#         if iter == max_iters
#             println("Warning: Maximum iterations reached:")
#             println("Final errors: convergence error = $convergence_error")
#             emissive     = updateSpectralEmission!(rtm, iter, F_matrices, sol_j,
#                                                    emitFrac, temperatures, emissive)
#             temperatures = updateTemperaturesSpectral!(rtm, emissive, emitFrac)
#         end
#     end

#     # ---- Write spectral results back to the domain ---------------------------
#     println("Writing spectral results to mesh...")
#     for bin in 1:K
#         bin_start = (bin - 1) * N + 1
#         bin_end   = bin * N
#         j_bin     = sol_j[bin_start:bin_end]

#         # e = D j ; r = j - e = R' j ; g_a = A' j = j - C j - r
#         e_bin   = (I - Diagonal(b[:,bin])*F_matrices[bin]') * j_bin
#         r_bin   = j_bin .- e_bin
#         g_a_bin = j_bin .- (I - F_matrices[bin]') * j_bin .- r_bin # C_matrices[bin]

#         writeResultsToDomain!(rtm, j_bin, g_a_bin, r_bin; spectral_bin=bin)
#     end

#     writeTemperaturesHeatSources!(rtm, temperatures)

#     # Per-bin energy conservation error (same field, same definition)
#     rtm.energy_error = G.([sum((I - F_matrices[i]') *
#                                 sol_j[(i - 1) * N + 1 : i * N]) for i in 1:K]) # C_matrices[i]

#     return rtm
# end