function buildSystemMatrices!(mesh::RayTracingDomain2D, F::Matrix{P}; 
                                spectral_bin::Int=1) where {P<:Real}
    
    # Count surfaces and volumes
    N_surfs = sum([sum(fine_face.solidWalls) for i in 1:length(mesh.coarse_mesh) 
                   for fine_face in mesh.fine_mesh[i]])
    N_vols = length([fine_face for i in 1:length(mesh.coarse_mesh) 
                     for fine_face in mesh.fine_mesh[i]])
    n = N_surfs + N_vols
    
    # Allocate output matrices
    B_save = zeros(P, n, n)
    S_infty_save = zeros(P, n, n)
    A_save = zeros(P, n, n)
    R_save = zeros(P, n, n)
    C_save = zeros(P, n, n)
    D_save = zeros(P, n, n)
    M_save = zeros(P, n, n)
    
    # Allocate workspace
    ws = TwoMatrixWorkspace{P}(N_surfs, N_vols)
    
    # Use work vectors for emissive powers (needed for M matrix)
    E_known = ws.work_vec1
    Q_known = ws.work_vec2
    
    # Step 1: Populate workspace from mesh
    populateWorkspace!(ws, mesh, spectral_bin)
    
    # Step 2: Compute emissive powers (needed for M matrix boundary conditions)
    computeEmissivePowersVariable!(mesh, ws, E_known, Q_known)
    
    # Step 3: Compute B matrix (reflectivity/scattering albedo)
    B = ws.matrix1
    fill!(B, zero(P))
    
    # Check if we need scattering/reflection calculations
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
    B_save .= B
    
    # Step 4: Compute K matrix (scattering kernel)
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
    
    # Step 5: Solve for S_infty (total exchange factor with multiple reflections)
    if !has_scattering
        copyto!(ws.matrix2, F)  # S_infty = F (no scattering)
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
        ws.matrix2 .= K \ F  # Solve (I - K) * S_infty = F
    end
    
    S_infty = ws.matrix2  # matrix2 now contains S_infty
    S_infty_save .= S_infty
    
    # Step 6: Compute A (absorption) and R (reflection) matrices
    for i in 1:n
        for j in 1:n
            # Get B values using helper function
            B_ji = getValueB(ws, j, N_surfs)
            B_ij = getValueB(ws, i, N_surfs)
            
            # A[i,j] = (1 - B[j,i]) * S_infty[i,j] * (1 - B[i,j])
            # R[i,j] = (1 - B[j,i]) * S_infty[i,j] * B[i,j]
            common_term = (1.0 - B_ji) * S_infty[i, j]
            A_save[i, j] = common_term * (1.0 - B_ij)
            R_save[i, j] = common_term * B_ij
        end
    end
    
    # Step 7: Compute C and D matrices
    C_save .= I - R_save' - A_save'
    D_save .= I - R_save'
    
    # Step 8: Build M matrix (system matrix for solving)
    M = ws.matrix1  # Reuse matrix1
    Q_vec_known = vcat(ws.bin_Qw_known, ws.bin_Qg_known)
    
    for i in 1:n
        for j in 1:n
            if Q_vec_known[i] == 1
                # Known heat flux boundary condition
                # M[i,j] = I[i,j] - A'[i,j] - R'[i,j] = I[i,j] - A[j,i] - R[j,i]
                if i == j
                    M[i, j] = 1.0 - A_save[i, j] - R_save[i, j]
                else
                    M[i, j] = -A_save[j, i] - R_save[j, i]
                end
            else
                # Known temperature boundary condition
                # M[i,j] = I[i,j] - R'[i,j] = I[i,j] - R[j,i]
                if i == j
                    M[i, j] = 1.0 - R_save[i, j]
                else
                    M[i, j] = -R_save[j, i]
                end
            end
        end
    end
    M_save .= M
    
    return C_save, D_save, M_save
end

function buildSystemMatrices!(domain::ViewFactorDomain3D, F::Matrix{P}; 
                                   spectral_bin::Int=1) where {P<:Real}
    
    # Count surfaces (no volumes in 3D!)
    N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    n = N_surfs
    
    # Allocate output matrices
    C_save = zeros(P, n, n)
    D_save = zeros(P, n, n)
    M_save = zeros(P, n, n)
    
    # Allocate workspace
    ws = SurfaceOnlyWorkspace{P}(N_surfs)
    
    # Use work vectors for emissive powers (needed for M matrix)
    E_known = ws.work_vec1
    Q_known = ws.work_vec2
    
    # Step 1: Populate workspace from domain
    populateWorkspace!(ws, domain, spectral_bin)
    
    # Step 2: Compute emissive powers (needed for M matrix boundary conditions)
    computeEmissivePowers!(E_known, Q_known, ws, N_surfs)
    
    # Step 3: Compute B matrix (surface reflectivity only!)
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
    
    # Step 4: Compute K matrix (reflection kernel)
    K = ws.matrix2
    for i in 1:n
        for j in 1:n
            K[i, j] = F[i, j] * B[i, j]
        end
    end
    
    # Step 5: Solve for S_infty (total exchange factor with multiple reflections)
    if !has_reflection
        copyto!(ws.matrix2, F)  # S_infty = F (no reflection)
    else
        # (I - K) * S_infty = F
        for i in 1:n
            for j in 1:n
                K[i, j] = (i == j ? 1.0 : 0.0) - K[i, j]
            end
        end
        ws.matrix2 .= K \ F  # Solve (I - K) * S_infty = F
    end
    
    # matrix2 now contains S_infty
    
    # Step 6: Compute A (absorption) and R (reflection) matrices
    A = zeros(P, n, n)
    R = zeros(P, n, n)
    for i in 1:n
        for j in 1:n
            # Get B values using helper function
            B_ji = getValueB(ws, j)
            B_ij = getValueB(ws, i)
            
            # A[i,j] = (1 - B[j,i]) * S_infty[i,j] * (1 - B[i,j])
            # R[i,j] = (1 - B[j,i]) * S_infty[i,j] * B[i,j]
            common_term = (1.0 - B_ji) * ws.matrix2[i, j]
            A[i, j] = common_term * (1.0 - B_ij)
            R[i, j] = common_term * B_ij
        end
    end
    
    # Step 7: Compute C and D matrices
    C_save .= I - R' - A'
    D_save .= I - R'
    
    # Step 8: Build M matrix (system matrix for solving)
    M_save = ws.matrix1  # Reuse matrix1
    Q_vec_known = ws.bin_Qw_known
    
    for i in 1:n
        for j in 1:n
            if Q_vec_known[i] == 1
                # Known heat flux boundary condition
                # M[i,j] = I[i,j] - A'[i,j] - R'[i,j] = I[i,j] - A[j,i] - R[j,i]
                if i == j
                    M_save[i, j] = 1.0 - A[i, j] - R[i, j]
                else
                    M_save[i, j] = -A[j, i] - R[j, i]
                end
            else
                # Known temperature boundary condition
                # M[i,j] = I[i,j] - R'[i,j] = I[i,j] - R[j,i]
                if i == j
                    M_save[i, j] = 1.0 - R[i, j]
                else
                    M_save[i, j] = -R[j, i]
                end
            end
        end
    end
    
    return C_save, D_save, M_save
end