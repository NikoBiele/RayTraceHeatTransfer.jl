
# Compute emissive powers (surfaces only!)
function computeEmissivePowers!(E_known::Vector{P}, Q_known::Vector{P}, 
                                    ws::SurfaceOnlyWorkspace{P}, N_surfs::Int) where {P}
    
    # Surface emissive powers
    for i in 1:N_surfs
        if ws.bin_Qw_known[i] == 0  # Known temperature
            E_known[i] = ws.epsw[i] * STEFAN_BOLTZMANN * ws.Area[i] * ws.Tw[i]^4
            Q_known[i] = 0
        else  # Known flux
            E_known[i] = 0
            Q_known[i] = ws.qw[i]
        end
    end
end

# Compute temperatures (surfaces only!)
function computeTemperatures!(T::Vector{P}, j::Vector{P}, r::Vector{P}, 
                                 ws::SurfaceOnlyWorkspace{P}, N_surfs::Int) where {P}
    
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
end

# Main grey solver for 3D surfaces - ALWAYS solves
function equilibriumSurfacesGrey3D!(domain::ViewFactorDomain3D, F::Matrix{P}; 
                                   spectral_bin::Int=1) where {P<:Real}
    println("=== 3D Surface-Only Grey Solver ===")
    
    # Count surfaces (no volumes!)
    N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    n = N_surfs
    
    println("Found $N_surfs surfaces")
    
    # Allocate workspace
    ws = SurfaceOnlyWorkspace{P}(N_surfs)
    
    # Use work vectors
    E_known = ws.work_vec1
    Q_known = ws.work_vec2
    
    # Step 1: Populate workspace
    println("Populating workspace...")
    populateWorkspace!(ws, domain, spectral_bin)
    
    # Step 2: Compute emissive powers
    println("Computing emissive powers...")
    computeEmissivePowers!(E_known, Q_known, ws, N_surfs)
    
    # Step 3: Compute B matrix (surface reflectivity only!)
    println("Computing B matrix...")
    B = ws.matrix1
    fill!(B, zero(P))
    
    # Check if we need reflection calculations
    has_reflection = sum(ws.epsw) < N_surfs
    
    if has_reflection
        # Surface reflectivity: B[i,j] = (1 - epsilon_j)
        for i in 1:n
            for j in 1:n
                B[i, j] = 1.0 - ws.epsw[j]
            end
        end
    end
    
    # Step 4: Compute K matrix
    println("Computing K matrix...")
    K = ws.matrix2
    for i in 1:n
        for j in 1:n
            K[i, j] = F[i, j] * B[i, j]
        end
    end
    
    # Step 5: Solve for S_infty
    println("Solving for S_infty...")
    
    if !has_reflection
        copyto!(ws.matrix2, F)  # S_infty = F
    else
        # (I - K) * S_infty = F
        for i in 1:n
            for j in 1:n
                K[i, j] = (i == j ? 1.0 : 0.0) - K[i, j]
            end
        end
        ws.matrix2 .= K \ F
    end
    
    # matrix2 now contains S_infty
    
    # Step 6: Build system matrix M
    println("Assembling linear system...")
    M = ws.matrix1  # Reuse matrix1
    b = ws.work_vec3
    
    Q_vec_known = ws.bin_Qw_known
    for i in 1:n
        b[i] = Q_vec_known[i] == 1 ? Q_known[i] : E_known[i]
    end
    
    # Build M using inline helper (clean and efficient!)
    for i in 1:n
        for j in 1:n
            B_ji = getValueB(ws, j)
            B_ij = getValueB(ws, i)
            
            common_term = (1.0 - B_ji) * ws.matrix2[i, j]
            A_ij = common_term * (1.0 - B_ij)
            R_ij = common_term * B_ij
            
            if Q_vec_known[i] == 1  # Known flux
                # Need A[j,i] and R[j,i]
                B_ij_swap = getValueB(ws, j)
                B_ji_swap = getValueB(ws, i)
                common_term_swap = (1.0 - B_ji_swap) * ws.matrix2[j, i]
                A_ji = common_term_swap * (1.0 - B_ij_swap)
                R_ji = common_term_swap * B_ij_swap
                
                M[i, j] = (i == j) ? (1.0 - A_ij - R_ij) : (-A_ji - R_ji)
            else  # Known temperature
                # Need R[j,i]
                B_ij_swap = getValueB(ws, j)
                B_ji_swap = getValueB(ws, i)
                common_term_swap = (1.0 - B_ji_swap) * ws.matrix2[j, i]
                R_ji = common_term_swap * B_ij_swap
                
                M[i, j] = (i == j) ? (1.0 - R_ij) : (-R_ji)
            end
        end
    end
    
    # Step 7: Solve
    println("Solving linear system...")
    j = M \ b
    
    # Step 8: Compute absorbed and reflected
    println("Computing absorbed and reflected energies...")
    Abs = ws.work_vec4
    r = ws.work_vec5
    fill!(Abs, zero(P))
    fill!(r, zero(P))
    
    for i in 1:n
        local_abs = zero(P)
        local_r = zero(P)
        
        B_i = getValueB(ws, i)
        
        for k in 1:n
            B_k = getValueB(ws, k)
            
            common_term = (1.0 - B_i) * ws.matrix2[k, i]
            A_ki = common_term * (1.0 - B_k)
            R_ki = common_term * B_k
            
            local_abs += A_ki * j[k]
            local_r += R_ki * j[k]
        end
        
        Abs[i] = local_abs
        r[i] = local_r
    end
    
    # Step 9: Compute temperatures
    println("Computing temperatures...")
    T = ws.work_vec1
    computeTemperatures!(T, j, r, ws, N_surfs)
    
    # Step 10: Write results to domain
    println("Writing results to domain...")
    writeResultsToDomain!(domain, j, Abs, r; T=T, spectral_bin=spectral_bin)

    # Step 11: Compute energy conservation error
    println("Computing energy conservation error...")
    domain.energy_error = sum(j - r - Abs)
    
    println("=== 3D Grey Solution Complete ===")
    
    return nothing
end

