# Fixed radiative transfer algorithm

# Updated compute_emissive_powers - use precomputed kappa
function computeEmissivePowersVariable!(mesh::RayTracingDomain2D, ws::TwoMatrixWorkspace{P},
                                        E_known::Vector{P}, Q_known::Vector{P}) where {P}
    
    # Build Q_vec_known
    Q_vec_known = vcat(ws.bin_Qw_known, ws.bin_Qg_known)
    
    # Surface emissive powers
    N_surfs = length(mesh.surface_mapping)
    for i in 1:N_surfs
        if Q_vec_known[i] == 0
            E_known[i] = ws.epsw[i] * STEFAN_BOLTZMANN * ws.Area[i] * ws.Tw[i]^4
            Q_known[i] = 0
        else
            E_known[i] = 0
            Q_known[i] = ws.qw[i]
        end
    end
    
    # Volume emissive powers - use precomputed kappa
    for vol_count in 1:length(mesh.volume_mapping)
        vol_idx = N_surfs + vol_count
        
        if Q_vec_known[vol_idx] == 0
            E_known[vol_idx] = 4 * ws.kappa_g[vol_count] * STEFAN_BOLTZMANN * ws.Volume[vol_count] * ws.Tg[vol_count]^4
            Q_known[vol_idx] = 0
        else
            E_known[vol_idx] = 0
            Q_known[vol_idx] = ws.qg[vol_count]
        end
    end
end

# Updated compute temperatures - use precomputed kappa
function computeTemperaturesVariable!(mesh::RayTracingDomain2D, ws::TwoMatrixWorkspace{P},
                                    T::Vector{P}, j::Vector{P}, r::Vector{P}) where {P}
    
    # Surface temperatures
    N_surfs = length(mesh.surface_mapping)
    for i in 1:N_surfs
        e_i = max(j[i] - r[i], 0.0)
        if ws.epsw[i] > 0.0 && ws.Area[i] > 0.0
            T[i] = (e_i / (ws.epsw[i] * STEFAN_BOLTZMANN * ws.Area[i]))^0.25
        else
            T[i] = 0.0
        end
        if isnan(T[i])
            T[i] = 0.0
        end
    end
    
    # Volume temperatures - use precomputed kappa
    for vol_count in 1:length(mesh.volume_mapping)
        vol_idx = N_surfs + vol_count
        
        e_i = max(j[vol_idx] - r[vol_idx], 0.0)
        if ws.kappa_g[vol_count] > 0.0 && ws.Volume[vol_count] > 0.0
            T[vol_idx] = (e_i / (4 * ws.kappa_g[vol_count] * ws.Volume[vol_count] * STEFAN_BOLTZMANN))^0.25
        else
            T[vol_idx] = 0.0
        end
        
        if isnan(T[vol_idx])
            T[vol_idx] = 0.0
        end
    end
end

# Main solver - UNCHANGED except using getValueB
function equilibriumGrey2D!(mesh::RayTracingDomain2D, F::Matrix{P}; spectral_bin::N=1) where {P<:Real, N<:Integer}
    println("=== Variable Extinction Memory-Optimized Steady State Solver ===")
    
    # Count surfaces and volumes
    N_surfs = length(mesh.surface_mapping)
    N_vols = length(mesh.volume_mapping)
    n = N_surfs + N_vols
        
    # Allocate workspace
    println("Allocating workspace...")
    ws = TwoMatrixWorkspace{P}(N_surfs, N_vols)
    
    # Use work vectors for main arrays
    E_known = ws.work_vec1
    Q_known = ws.work_vec2
    
    # Step 1: Populate workspace from mesh
    println("Populating workspace from mesh...")
    populateWorkspace!(ws, mesh, spectral_bin)
    
    # Step 2: Compute emissive powers with variable extinction
    println("Computing emissive powers with variable extinction...")
    computeEmissivePowersVariable!(mesh, ws, E_known, Q_known)
    
    # Step 3: Compute B matrix with variable extinction
    println("Computing B matrix with variable extinction...")
    B = ws.matrix1
    fill!(B, zero(P))
    
    # Check if we need scattering calculations
    has_scattering = any(ws.omega_g .> 1e-6)
    has_reflection = sum(ws.epsw) < n
    
    if has_scattering || has_reflection
        # Surface reflectivity terms
        for i in 1:n
            for j in 1:N_surfs
                B[i, j] = 1.0 - ws.epsw[j]
            end
        end
        
        # Volume scattering terms - use precomputed omega
        for vol_count in 1:N_vols
            vol_idx = N_surfs + vol_count
            for i in 1:n
                B[i, vol_idx] = ws.omega_g[vol_count]
            end
        end
    end
    
    # Step 4: Compute K matrix in matrix2
    println("Computing K matrix...")
    K = ws.matrix2
    for i in 1:n
        for j in 1:n
            if B[i, j] != 0.0 && F[i, j] != 0.0
                K[i, j] = F[i, j] * B[i, j]
            else
                K[i, j] = 0.0
            end
        end
    end
    
    # Step 5: Solve for S_infty (overwrite K in matrix2)
    println("Solving for S_infty...")
    
    if !has_scattering
        copyto!(ws.matrix2, F)  # S_infty = F
    else
        # Convert K to (I - K) in matrix2
        for i in 1:n
            for j in 1:n
                if i == j
                    K[i, j] = 1.0 - K[i, j]
                else
                    K[i, j] = -K[i, j]
                end
            end
        end
        ws.matrix2 .= K \ F  # Solve and store in matrix2
    end
    
    # matrix2 now contains S_infty
    
    # Step 6: Build system matrix M directly in matrix1 (overwrite B)
    println("Assembling linear system...")
    M = ws.matrix1  # Reuse matrix1 for system matrix
    b = ws.work_vec3
    
    # Build Q_vec_known and RHS
    Q_vec_known = vcat(ws.bin_Qw_known, ws.bin_Qg_known)
    for i in 1:n
        if Q_vec_known[i] == 1
            b[i] = Q_known[i]
        else
            b[i] = E_known[i]
        end
    end
    
    # Build system matrix M - NOW using getValueB instead of get_local_omega
    for i in 1:n
        for j in 1:n
            # Get B values using helper (no nested loops!)
            B_ji = getValueB(ws, j, N_surfs)
            B_ij = getValueB(ws, i, N_surfs)
            
            # A[i,j] = (1 - B[j,i]) * S_infty[i,j] * (1 - B[i,j])
            # R[i,j] = (1 - B[j,i]) * S_infty[i,j] * B[i,j]
            common_term = (1.0 - B_ji) * ws.matrix2[i, j]
            A_ij = common_term * (1.0 - B_ij)
            R_ij = common_term * B_ij
            
            # Build system matrix based on boundary condition
            if Q_vec_known[i] == 1
                # M[i,j] = I[i,j] - A'[i,j] - R'[i,j] = I[i,j] - A[j,i] - R[j,i]
                if i == j
                    M[i, j] = 1.0 - A_ij - R_ij
                else
                    # Need A[j,i] and R[j,i] - compute with swapped indices
                    B_ij_swap = getValueB(ws, j, N_surfs)
                    B_ji_swap = getValueB(ws, i, N_surfs)
                    common_term_swap = (1.0 - B_ji_swap) * ws.matrix2[j, i]
                    A_ji = common_term_swap * (1.0 - B_ij_swap)
                    R_ji = common_term_swap * B_ij_swap
                    M[i, j] = -A_ji - R_ji
                end
            else
                # M[i,j] = I[i,j] - R'[i,j] = I[i,j] - R[j,i]
                if i == j
                    M[i, j] = 1.0 - R_ij
                else
                    # Need R[j,i]
                    B_ij_swap = getValueB(ws, j, N_surfs)
                    B_ji_swap = getValueB(ws, i, N_surfs)
                    common_term_swap = (1.0 - B_ji_swap) * ws.matrix2[j, i]
                    R_ji = common_term_swap * B_ij_swap
                    M[i, j] = -R_ji
                end
            end
        end
    end
    
    # Step 7: Solve linear system
    println("Solving linear system...")
    j = M \ b
    
    # Step 8: Compute Abs = A' * j and r = R' * j
    println("Computing absorbed and reflected energies...")
    Abs = ws.work_vec4
    r = ws.work_vec5
    fill!(Abs, zero(P))
    fill!(r, zero(P))
    
    # Compute matrix-vector products - NOW using getValueB
    for i in 1:n
        local_abs = zero(P)
        local_r = zero(P)
        
        for k in 1:n
            # Get B values using helper
            B_ik = getValueB(ws, i, N_surfs)
            B_ki = getValueB(ws, k, N_surfs)
            
            common_term = (1.0 - B_ik) * ws.matrix2[k, i]
            A_ki = common_term * (1.0 - B_ki)
            R_ki = common_term * B_ki
            
            local_abs += A_ki * j[k]
            local_r += R_ki * j[k]
        end
        
        Abs[i] = local_abs
        r[i] = local_r
    end
    
    # Step 9: Compute temperatures with variable extinction
    println("Computing temperatures with variable extinction...")
    T = ws.work_vec1  # Reuse work_vec1
    computeTemperaturesVariable!(mesh, ws, T, j, r)
    
    # Step 10: Write results to mesh
    println("Writing results to mesh...")
    writeResultsToDomain!(mesh, j, Abs, r; T=T, spectral_bin=spectral_bin)

    # Step 11: Compute energy conservation error
    println("Computing energy conservation error...")
    mesh.energy_error = sum(j - r - Abs)
    
    println("=== Variable Extinction Steady State Solution Complete ===")
    
    return nothing
end