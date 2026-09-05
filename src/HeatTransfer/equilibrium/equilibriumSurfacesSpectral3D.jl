function equilibriumSurfacesSpectral3D!(domain::SurfaceDomain3D{G,P}, F::AbstractMatrix; 
                               max_iters::Int=1000, convergence_tol::G=1e-12, verbose::Bool=true) where
                               {G,P<:Integer}
    """
    Solve 3D spectral surface radiation using the optimal solver.
    
    Automatically chooses between:
    - Direct emission solver: for ε=1 everywhere (orders of magnitude faster for many bins)
    - Full iterative solver: for general problems with reflection
    
    Set force_full=true to disable optimization and always use full solver.
    """
    
    # Check if we can use the fast direct solver
    if domain.uniform_epsilon
        b_check = get_b(domain)
        if any(b_check .> 1e-12)
            error("""
            uniform_epsilon dispatch selected the direct solver, but the
            domain has reflecting surfaces (max b = $(maximum(b_check))).
            The direct solver is not valid for this configuration.
            This indicates uniform_epsilon was set inconsistently with the
            surface properties; the automatic detection selects the correct
            solver - avoid setting it manually.
            """)
        end
        verbose && println("=== Using DIRECT solver ===")
        return equilibriumSurfacesSpectral3D_direct!(domain, F; verbose=verbose)
    else
        verbose && println("=== Using FULL solver ===")
        return equilibriumSurfacesSpectral3D_woodbury!(domain, F; max_iters=max_iters,
                                                convergence_tol=convergence_tol, verbose=verbose)
    end
end

function equilibriumSurfacesSpectral3D_woodbury!(domain::SurfaceDomain3D{G,P}, F::AbstractMatrix;
                               max_iters::Int=1000,
                               convergence_tol::G=1e-12, verbose::Bool=true) where {G,P<:Integer}

    if isnothing(domain.wavelength_band_limits)
        error("""
        Spectral solve requires wavelength band limits to be set.
 
        Please set mesh.wavelength_band_limits before calling steadyStateSpectral3D!:
 
        Example, using logarithmic spacing:
            mesh.wavelength_band_limits = 10 .^ range(log10(0.0000001), log10(0.001), length=51)
        """)
    end
    if length(domain.wavelength_band_limits) < 4
        error("wavelength_band_limits must have at least 4 values (defining 3 bins)")
    end
    if any(domain.wavelength_band_limits .<= 0)
        error("wavelength_band_limits must all be positive (wavelengths > 0)")
    end
    if any(diff(domain.wavelength_band_limits) .<= 0)
        error("wavelength_band_limits must be strictly increasing (no duplicates)")
    end
 
    verbose && println("=== 3D Spectral Surface Radiation Solver (WOODBURY) ===")
    verbose && println("Spectral mode: $(domain.spectral_mode)")
    verbose && println("Number of spectral bins: $(domain.n_spectral_bins)")
 
    K = domain.n_spectral_bins
    N = sum([length(superface.subFaces) for superface in domain.facesMesh])
 
    verbose && println("\nComputing GERT matrices for each spectral band...")
    verbose && println("(Using same view factor matrix F for all bands)")
 
    b = get_b(domain)
    b_max = maximum(b)
    if b_max > 1 - 1e-6
        @warn """
        b_max = $b_max approaches 1 (near-perfectly reflecting enclosure).
        The unridged Cholesky of D'D is certified for b_max < 1 with
        conditioning ~ (1 - b_max)^-2; this configuration is at the edge
        of that certificate and accuracy may degrade.
        """
    end
 
    Ft = Matrix(F')                       # dense throughout (3D view factor F is dense)
    M_matrices = Vector{Matrix{G}}(undef, K)
    D_matrices = Vector{Matrix{G}}(undef, K)
    for k in 1:K
        M_matrices[k] = Matrix(buildSystemMatrix(domain, F; spectral_bin=k))
        D_matrices[k] = Matrix{G}(I, N, N) - Diagonal(b[:, k]) * Ft
    end
 
    # ---- Woodbury precomputation: per-bin Cholesky + capacitance -------------
    verbose && println("Factorizing per-bin normal matrices + capacitance...")
    DtD_factors = Vector{Cholesky{G, Matrix{G}}}(undef, K)
    DtDinv_Mt   = Vector{Matrix{G}}(undef, K)
    inner = Matrix{G}(I, N, N)
    for k in 1:K
        Dk   = D_matrices[k]
        chol = cholesky(Symmetric(Dk' * Dk))          # no ridge: certified SPD
        DtD_factors[k] = chol
        DtDinv_Mt[k]   = chol \ Matrix(M_matrices[k]')
        inner .+= M_matrices[k] * DtDinv_Mt[k]
    end
    inner_factor = cholesky(Symmetric(inner))
 
    # ---- Boundary conditions and iteration state -----------------------------
    boundary, temperatures, emissive = setupBoundaryConditions(domain)
 
    sol_j          = zeros(G, K * N)
    previous_sol_j = zeros(G, K * N)
    emitFrac       = getWeightedEmissionFractions(domain, temperatures)
 
    u_list   = [zeros(G, N) for _ in 1:K]
    Mu_total = zeros(G, N)
 
    nthr    = Threads.nthreads()
    chunks  = collect(Iterators.partition(1:K, cld(K, nthr)))
    nchunks = length(chunks)
    Mu_partial = [zeros(G, N) for _ in 1:nchunks]
    e_k_buf    = [zeros(G, N) for _ in 1:nchunks]
    rhs_k_buf  = [zeros(G, N) for _ in 1:nchunks]
 
    verbose && println("\nStarting spectral steady-state iteration (Woodbury)...")
    for iter = 1:max_iters
 
        # ---- Woodbury solve of the bordered LS normal equations -----------
        for buf in Mu_partial; fill!(buf, zero(G)); end
        Threads.@sync for (t, chunk) in enumerate(chunks)
            Threads.@spawn begin
                e_k = e_k_buf[t]; rhs_k = rhs_k_buf[t]; Mu_t = Mu_partial[t]
                for k in chunk
                    Mk = M_matrices[k]; Dk = D_matrices[k]
                    @views e_k .= emissive .* emitFrac[:, k]
                    mul!(rhs_k, Mk', boundary)
                    mul!(rhs_k, Dk', e_k, one(G), one(G))
                    ldiv!(u_list[k], DtD_factors[k], rhs_k)
                    mul!(Mu_t, Mk, u_list[k], one(G), one(G))
                end
            end
        end
        fill!(Mu_total, zero(G))
        for t in 1:nchunks; Mu_total .+= Mu_partial[t]; end
        w = inner_factor \ Mu_total
        Threads.@sync for (t, chunk) in enumerate(chunks)
            Threads.@spawn begin
                corr = rhs_k_buf[t]
                for k in chunk
                    mul!(corr, DtDinv_Mt[k], w)
                    @views sol_j[(k-1)*N+1 : k*N] .= u_list[k] .- corr
                end
            end
        end
 
        # ---- emission / temperature update (same order as 2D Woodbury) ----
        emissive = updateSpectralEmission!(domain, iter, F, sol_j,
                                           emitFrac, temperatures, emissive)
        temperatures = updateTemperaturesSpectral!(domain, emissive,
                                           getBinsEmissionFractions(domain, temperatures))
        emitFrac = getWeightedEmissionFractions(domain, temperatures)
 
        # ---- convergence check (same criterion as _full!) ------------------
        convergence_error = norm(sol_j - previous_sol_j) / norm(sol_j)
        previous_sol_j .= sol_j
        if iter % 20 == 0
            verbose && println("Iteration $iter: convergence error = $convergence_error")
        end
 
        if iter > 1 && convergence_error < convergence_tol
            verbose && println("Converged after $iter iterations")
            domain.energy_error = G.([sum((I - F') * sol_j[(i - 1) * N + 1:i * N])
                               for i in 1:K])
            verbose && println("Energy conservation errors by band: $(domain.energy_error)")
            break
        end
 
        if iter == max_iters
            @warn "Warning: Maximum iterations reached. Final error = $convergence_error"
        end
    end
 
    # ---- Write results to mesh (verbatim from _full!) ------------------------
    verbose && println("\nWriting spectral results to mesh...")
    for bin in 1:K
        bin_start = (bin - 1) * N + 1
        bin_end = bin * N
        j_bin = sol_j[bin_start:bin_end]
 
        e_bin = (I - Diagonal(b[:,bin]) * F') * j_bin     # e = D * j
        r_bin = j_bin - e_bin                             # r = j - e = R' * j
        g_a_bin = j_bin - (I - F') * j_bin - r_bin        # g_a = A' * j
 
        surf_count = 0
        for superface in domain.facesMesh
            for subface in superface.subFaces
                surf_count += 1
                if isa(subface.e_w, Vector)
                    subface.j_w[bin] = j_bin[surf_count]
                    subface.e_w[bin] = max(e_bin[surf_count], 0.0)
                    subface.r_w[bin] = max(r_bin[surf_count], 0.0)
                    subface.g_a_w[bin] = max(g_a_bin[surf_count], 0.0)
                    subface.g_w[bin] = subface.g_a_w[bin] + subface.r_w[bin]
                    subface.i_w[bin] = j_bin[surf_count] / (π * subface.area)
                end
            end
        end
    end
 
    writeTemperaturesHeatSources!(domain, temperatures)
 
    verbose && println("=== 3D Spectral Solution Complete (WOODBURY) ===")
end

function equilibriumSurfacesSpectral3D_direct!(domain::SurfaceDomain3D{G,P}, F::AbstractMatrix;
                                                verbose::Bool=true) where {G,P<:Integer}
    """
    OPTIMIZED direct emission solver for 3D non-scattering, non-reflecting problems.
    
    Key optimization: Instead of solving (n_bins+1)*N × n_bins*N block system for j,
    we solve N×N system for e directly, then recover j from emission fractions.
    
    Requires: ε=1 everywhere (no reflection) - 3D has no volumes so no scattering check needed
    """
    
    # Validate spectral setup at entry
    if isnothing(domain.wavelength_band_limits)
        error("""
        Spectral solve requires wavelength band limits to be set.
        
        Please set mesh.wavelength_band_limits before calling steadyStateSpectral3D_direct!:
        
        Example, using logarithmic spacing:
            mesh.wavelength_band_limits = 10 .^ range(log10(0.0000001), log10(0.001), length=51)
        """)
    end
    
    # Additional validation
    if length(domain.wavelength_band_limits) < 4
        error("wavelength_band_limits must have at least 4 values (defining 3 bins)")
    end
    
    if any(domain.wavelength_band_limits .<= 0)
        error("wavelength_band_limits must all be positive (wavelengths > 0)")
    end

    if any(diff(domain.wavelength_band_limits) .<= 0)
        error("wavelength_band_limits must be strictly increasing (no duplicates)")
    end

    verbose && println("=== 3D Spectral Surface Radiation Solver (DIRECT) ===")
    verbose && println("Spectral mode: $(domain.spectral_mode)")
    verbose && println("Number of spectral bins: $(domain.n_spectral_bins)")
    verbose && println("Using optimized direct emission solver")
    
    # Get system size
    N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    
    verbose && println("\nComputing GERT matrices for each spectral band...")
    verbose && println("(Using same view factor matrix F for all bands)")
    
    M = buildSystemMatrix(domain, F[1:N_surfs,1:N_surfs]; spectral_bin=1) # C, D
    
    # Setup boundary conditions
    verbose && println("Setting up boundary conditions...")
    boundary, temperatures, emissive = setupBoundaryConditions(domain)
    
    # Iteration loop
    verbose && println("\nStarting spectral direct solve...")
    # Build reduced system: A_reduced = sum_bin M_bin * Diagonal(emitFrac[:, bin])
    # This is the key optimization! N×N system instead of (n_bins*N)×(n_bins*N)        
    # DIRECT solve
    j_tot = M \ boundary
        
    # Update emissive powers and temperatures
    b = get_b(domain)
    emissive = (I - Diagonal(b[:,1])*F[1:N_surfs,1:N_surfs]')*j_tot

    # iterate to find the correct blackbody spectral fractions
    emitFrac = getBinsEmissionFractions(domain, temperatures)
    temperatures = updateTemperaturesSpectral!(domain, emissive, emitFrac)
    temperatures_prev = temperatures
    for iter = 1:10
        emitFrac = getBinsEmissionFractions(domain, temperatures)
        temperatures = updateTemperaturesSpectral!(domain, emissive, emitFrac)
        max_iteration_change = maximum(abs.(temperatures - temperatures_prev))
        if max_iteration_change < 1e-3
            break
        end
        temperatures_prev = temperatures
    end
                    
    # Recover j for each bin: j_bin = emitFrac[:, bin] .* e
    sol_j = zeros(G, domain.n_spectral_bins * N_surfs)
    weightedFrac = getWeightedEmissionFractions(domain, temperatures)
    for bin in 1:domain.n_spectral_bins
        bin_start = (bin - 1) * N_surfs + 1
        bin_end = bin * N_surfs
        sol_j[bin_start:bin_end] = weightedFrac[:, bin] .* j_tot
    end

    # Compute energy conservation error for each band
    domain.energy_error = G.([sum((I - F[1:N_surfs,1:N_surfs]') * sol_j[(i - 1) * N_surfs + 1:i * N_surfs]) 
                        for i in 1:domain.n_spectral_bins]) # C
    verbose && println("Energy conservation errors by band: $(domain.energy_error)")

    for bin in 1:domain.n_spectral_bins
        # Extract solution for this band
        bin_start = (bin - 1) * N_surfs + 1
        bin_end = bin * N_surfs
        j_bin = sol_j[bin_start:bin_end]
        
        # Compute e, r, g_a from GERT matrices and radiosity
        e_bin = (I - Diagonal(b[:,1])*F[1:N_surfs,1:N_surfs]') * j_bin           # e = D * j
        r_bin = j_bin - e_bin                     # r = j - e = R' * j
        g_a_bin = j_bin - (I - F[1:N_surfs,1:N_surfs]') * j_bin - r_bin  # g_a = A' * j = j - C*j - r # C
        
        # Write to mesh for this band
        surf_count = 0
        for superface in domain.facesMesh
            for subface in superface.subFaces
                surf_count += 1
                
                if isa(subface.e_w, Vector)
                    subface.j_w[bin] = j_bin[surf_count]
                    subface.e_w[bin] = max(e_bin[surf_count], 0.0)
                    subface.r_w[bin] = max(r_bin[surf_count], 0.0)
                    subface.g_a_w[bin] = max(g_a_bin[surf_count], 0.0)
                    subface.g_w[bin] = subface.g_a_w[bin] + subface.r_w[bin]
                    subface.i_w[bin] = j_bin[surf_count] / (π * subface.area)
                end
            end
        end
    end

    # Write temperatures and heat sources
    writeTemperaturesHeatSources!(domain, temperatures)
    
    verbose && println("=== 3D Spectral Solution Complete (DIRECT) ===")
end

###### ORIGINAL SPECTRAL SOLVER - kept for reference ######

# function equilibriumSurfacesSpectral3D_full!(domain::ViewFactorDomain3D{G,P}, F::AbstractMatrix; 
#                                max_iters::Int=1000,
#                                convergence_tol::G=1e-12, verbose::Bool=true) where {G,P<:Integer}
    
#     # Validate spectral setup at entry
#     if isnothing(domain.wavelength_band_limits)
#         error("""
#         Spectral solve requires wavelength band limits to be set.
        
#         Please set mesh.wavelength_band_limits before calling steadyStateSpectral3D!:
        
#         Example, using logarithmic spacing:
#             mesh.wavelength_band_limits = 10 .^ range(log10(0.0000001), log10(0.001), length=51)
#         """)
#     end
    
#     # Additional validation
#     if length(domain.wavelength_band_limits) < 4
#         error("wavelength_band_limits must have at least 4 values (defining 3 bins)")
#     end
    
#     if any(domain.wavelength_band_limits .<= 0)
#         error("wavelength_band_limits must all be positive (wavelengths > 0)")
#     end

#     if any(diff(domain.wavelength_band_limits) .<= 0)
#         error("wavelength_band_limits must be strictly increasing (no duplicates)")
#     end

#     verbose && println("=== 3D Spectral Surface Radiation Solver ===")
#     verbose && println("Spectral mode: $(domain.spectral_mode)")
#     verbose && println("Number of spectral bins: $(domain.n_spectral_bins)")
    
#     # Get system size
#     N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    
#     verbose && println("\nComputing GERT matrices for each spectral band...")
#     verbose && println("(Using same view factor matrix F for all bands)")
    
#     # Build system matrices for each spectral band
#     # Each band has different epsilon, so different M
#     M_matrices = Vector{AbstractMatrix}()
    
#     for bin in 1:domain.n_spectral_bins
#         verbose && println("  Building matrices for spectral bin $bin...")
#         M = buildSystemMatrix(domain, F; spectral_bin=bin)
#         push!(M_matrices, M)
#     end
    
#     verbose && println("\nAssembling block matrix structure...")
#     # Build block matrix (same structure as 2D)
#     if F isa SparseMatrixCSC
#         block_matrix = spzeros((domain.n_spectral_bins + 1) * N_surfs, 
#                         domain.n_spectral_bins * N_surfs)
#     else
#         block_matrix = zeros(G, (domain.n_spectral_bins + 1) * N_surfs, 
#                             domain.n_spectral_bins * N_surfs)
#     end
#     b = get_b(domain)
    
#     for i = 1:(domain.n_spectral_bins + 1)
#         for j = 1:domain.n_spectral_bins
#             row_start = (i - 1) * N_surfs + 1
#             row_end = i * N_surfs
#             col_start = (j - 1) * N_surfs + 1
#             col_end = j * N_surfs
            
#             if i == 1
#                 block_matrix[row_start:row_end, col_start:col_end] = M_matrices[j]
#             elseif i == j + 1
#                 block_matrix[row_start:row_end, col_start:col_end] = I - Diagonal(b[:,j])*F'
#             end
#         end
#     end
#     # Outside iteration loop - factorize once
#     Factorization = qr(block_matrix)
    
#     # Setup boundary conditions
#     verbose && println("Setting up boundary conditions...")
#     boundary, temperatures, emissive = setupBoundaryConditions(domain)
    
#     # convergence tolerance
#     convergenceTolerance = G == BigFloat ? sqrt(eps(G)) : 1e-11

#     # Iteration loop
#     sol_j = zeros(G, domain.n_spectral_bins * N_surfs)
#     previous_sol_j = zeros(G, domain.n_spectral_bins * N_surfs)
#     emitFrac = getBinsEmissionFractions(domain, temperatures)
        
#     verbose && println("\nStarting spectral iteration...")
#     for iter = 1:max_iters
        
#         # Update emissive powers and temperatures
#         emissive = updateSpectralEmission!(domain, iter, F, sol_j, 
#                                              emitFrac, temperatures, emissive)
#         temperatures = updateTemperaturesSpectral!(domain, emissive,
#                                              getBinsEmissionFractions(domain, temperatures))
#         emitFrac = getWeightedEmissionFractions(domain, temperatures)
        
#         # Solve block system
#         emissive_pow_vec = reduce(vcat, [emissive for _ in 1:domain.n_spectral_bins])
#         b_e_matrix = Diagonal(reduce(vcat, [boundary; emissive_pow_vec]))
#         sol_j .= Factorization \ (b_e_matrix*[ones(G, N_surfs); emitFrac[:]])
        
#         # Check convergence
#         convergence_error = norm(sol_j - previous_sol_j) / norm(sol_j)
#         previous_sol_j .= sol_j
#         if iter % 20 == 0
#             verbose && println("Iteration $iter: convergence error = $convergence_error")
#         end
        
#         verbose && println("Iteration $iter: convergence error = $convergence_error")
        
#         if iter > 1 && convergence_error < convergence_tol
#             verbose && println("Converged after $iter iterations")
            
#             # Compute energy conservation error for each band
#             domain.energy_error = G.([sum((I - F') * sol_j[(i - 1) * N_surfs + 1:i * N_surfs]) 
#                                for i in 1:domain.n_spectral_bins]) # C_matrices[i]
#             verbose && println("Energy conservation errors by band: $(domain.energy_error)")
            
#             break
#         end
        
#         if iter == max_iters
#             @warn "Warning: Maximum iterations reached. Final error = $convergence_error"
#         end
#     end
    
#     verbose && println("\nWriting spectral results to mesh...")
#     for bin in 1:domain.n_spectral_bins
#         # Extract solution for this band
#         bin_start = (bin - 1) * N_surfs + 1
#         bin_end = bin * N_surfs
#         j_bin = sol_j[bin_start:bin_end]
        
#         # Compute e, r, g_a from GERT matrices and radiosity
#         e_bin = (I - Diagonal(b[:,bin])*F') * j_bin           # e = D * j
#         r_bin = j_bin - e_bin                     # r = j - e = R' * j
#         g_a_bin = j_bin - (I - F') * j_bin - r_bin  # g_a = A' * j = j - C*j - r # C_matrices[bin]
        
#         # Write to mesh for this band
#         surf_count = 0
#         for superface in domain.facesMesh
#             for subface in superface.subFaces
#                 surf_count += 1
                
#                 if isa(subface.e_w, Vector)
#                     subface.j_w[bin] = j_bin[surf_count]
#                     subface.e_w[bin] = max(e_bin[surf_count], 0.0)
#                     subface.r_w[bin] = max(r_bin[surf_count], 0.0)
#                     subface.g_a_w[bin] = max(g_a_bin[surf_count], 0.0)
#                     subface.g_w[bin] = subface.g_a_w[bin] + subface.r_w[bin]
#                     subface.i_w[bin] = j_bin[surf_count] / (π * subface.area)
#                 end
#             end
#         end
#     end
    
#     # Write temperatures and heat sources
#     writeTemperaturesHeatSources!(domain, temperatures)
    
#     verbose && println("=== 3D Spectral Solution Complete ===")
# end