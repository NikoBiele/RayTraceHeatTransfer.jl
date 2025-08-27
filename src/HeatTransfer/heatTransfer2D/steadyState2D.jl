function steadyState_simple!(mesh::RayTracingMeshOptim, F::Matrix{P}, gas::GasProperties) where {P<:Real}
    # Count total number of surfaces and volumes
    N_surfs = sum([sum(fine_face.solidWalls) for i in 1:length(mesh.coarse_mesh) for fine_face in mesh.fine_mesh[i]])
    N_vols = length([fine_face for i in 1:length(mesh.coarse_mesh) for fine_face in mesh.fine_mesh[i]])
    println("Found $N_surfs surfaces and $N_vols volumes.")

    println("Initializing sparse arrays...")
    beta_sparse_lim = 2.0
    Area = gas.beta > beta_sparse_lim ? spzeros(P, N_surfs) : zeros(P, N_surfs)
    epsw = gas.beta > beta_sparse_lim ? spzeros(P, N_surfs) : zeros(P, N_surfs)
    Tw = gas.beta > beta_sparse_lim ? spzeros(P, N_surfs) : zeros(P, N_surfs)
    qw = gas.beta > beta_sparse_lim ? spzeros(P, N_surfs) : zeros(P, N_surfs)
    bin_Qw_known = gas.beta > beta_sparse_lim ? spzeros(Int, N_surfs) : zeros(Int, N_surfs)
    Volume = gas.beta > beta_sparse_lim ? spzeros(P, N_vols) : zeros(P, N_vols)
    Tg = gas.beta > beta_sparse_lim ? spzeros(P, N_vols) : zeros(P, N_vols)
    qg = gas.beta > beta_sparse_lim ? spzeros(P, N_vols) : zeros(P, N_vols)
    bin_Qg_known = gas.beta > beta_sparse_lim ? spzeros(Int, N_vols) : zeros(Int, N_vols)

    surf_count = 0
    vol_count = 0
    println("Populating arrays...")
    for i in 1:length(mesh.coarse_mesh)
        for (k, fine_face) in enumerate(mesh.fine_mesh[i])
            vol_count += 1
            Volume[vol_count] = fine_face.volume
            Tg[vol_count] = fine_face.T_in_g
            qg[vol_count] = fine_face.q_in_g
            bin_Qg_known[vol_count] = fine_face.T_in_g < 0.0 ? 1 : 0

            for (wall_idx, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surf_count += 1
                    Area[surf_count] = fine_face.area[wall_idx]
                    epsw[surf_count] = fine_face.epsilon[wall_idx]
                    Tw[surf_count] = fine_face.T_in_w[wall_idx]
                    qw[surf_count] = fine_face.q_in_w[wall_idx]
                    bin_Qw_known[surf_count] = fine_face.T_in_w[wall_idx] < 0.0 ? 1 : 0
                end
            end
        end
    end

    println("Calculating sources and emissive powers...")
    Q_vec_known = vcat(bin_Qw_known, bin_Qg_known)

    sigma = 5.670374419e-8  # Stefan-Boltzmann constant
    Ew_known = gas.beta > beta_sparse_lim ? spzeros(P, N_surfs) : zeros(P, N_surfs)
    Qw_known = gas.beta > beta_sparse_lim ? spzeros(P, N_surfs) : zeros(P, N_surfs)
    for i in 1:N_surfs
        if Q_vec_known[i] == 0
            Ew_known[i] = epsw[i] * sigma * Area[i] * Tw[i]^4
        else
            Qw_known[i] = qw[i] * Area[i]
        end
    end

    Eg_known = gas.beta > beta_sparse_lim ? spzeros(P, N_vols) : zeros(P, N_vols)
    Qg_known = gas.beta > beta_sparse_lim ? spzeros(P, N_vols) : zeros(P, N_vols)
    for i in 1:N_vols
        if Q_vec_known[i+N_surfs] == 0
            Eg_known[i] = 4 * gas.kappa * sigma * Volume[i] * Tg[i]^4
        else
            Qg_known[i] = qg[i] * Volume[i]
        end
    end

    E_known = vcat(Ew_known, Eg_known)
    Q_known = vcat(Qw_known, Qg_known)

    function compute_AR(F, epsw, omega, N_surfs, N_vols)
        n = N_surfs + N_vols
        B = gas.beta > beta_sparse_lim ? spzeros(P, n, n) : zeros(P, n, n)
        K = gas.beta > beta_sparse_lim ? spzeros(P, n, n) : zeros(P, n, n)

        if gas.omega > 1e-6 || sum(epsw) < n
            for i = 1:n
                for j = 1:N_surfs
                    B[i,j] = 1.0-epsw[j]
                end
                for j = N_surfs+1:n
                    B[i,j] = omega
                end
            end
            for i = 1:n
                for j = 1:n
                    # Only compute and store if result will be non-zero
                    if B[i,j] != 0.0 && F[i,j] != 0.0
                        K[i,j] = F[i,j]*B[i,j]
                    end
                end
            end
        end
        
        println("Solving for S_infty...")
        S_infty = isapprox(gas.sigma_s, 0.0, atol=1e-10) ? F : (I-K)\F

        println("Building A and R matrices...")
        A = gas.beta > beta_sparse_lim ? spzeros(P, n, n) : zeros(P, n, n)
        R = gas.beta > beta_sparse_lim ? spzeros(P, n, n) : zeros(P, n, n)
        if gas.beta > beta_sparse_lim
            for i in 1:n
                for j in 1:n
                    if S_infty[i,j] > 1e-9 && B[i] < 0.999 && B[j] < 0.999
                        A[i,j] = (1 - B[j,i]) * S_infty[i,j] * (1 - B[i,j])
                        R[i,j] = (1 - B[j,i]) * S_infty[i,j] * B[i,j]
                    end
                end
            end
        else
            A = (1 .- B)' .* S_infty .* (1 .- B)
            R = (1 .- B)' .* S_infty .* B
        end

        return B, S_infty, A, R
    end

    println("Computing A and R matrices...")
    B, S_infty, A, R = compute_AR(F, epsw, gas.omega, N_surfs, N_vols)

    println("Assembling linear system...")
    n = N_surfs + N_vols
    C = I - A' - R'
    D = I - R'
    M = gas.beta > beta_sparse_lim ? spzeros(P, n, n) : zeros(P, n, n)
    b = zeros(P, n)
    for i in 1:n
        if Q_vec_known[i] == 1
            b[i] = Q_known[i]
            M[i,:] = C[i,:]
        else
            b[i] = E_known[i]
            M[i,:] = D[i,:]
        end
    end

    println("Solving for j...")
    j = M \ b

    println("Calculating absorbed energy, emitted energy, and other variables...")
    Abs = A' * j
    r = R' * j
    g = Abs + r
    e = max.(j .- r, 0.0)
    println("sum Q = 1^T*C*j = ", sum(C*j) )
    Q = e - Abs

    T = gas.beta > beta_sparse_lim ? spzeros(P, n) : zeros(P, n)
    for i in 1:n
        if i <= N_surfs
            T[i] = (e[i] / (epsw[i] * sigma * Area[i]))^0.25
        else
            T[i] = (e[i] / (4 * gas.kappa * Volume[i-N_surfs] * sigma))^0.25
        end
    end
    T[occursin.("NaN", string.(T))] .= 0.0

    println("Writing results back into the mesh...")
    surf_count = 0
    vol_count = 0
    for (i, coarse_face) in enumerate(mesh.coarse_mesh)
        for (k, fine_face) in enumerate(mesh.fine_mesh[i])
            vol_count += 1
            fine_face.T_g = T[N_surfs + vol_count]
            fine_face.j_g = j[N_surfs + vol_count]
            fine_face.g_a_g = Abs[N_surfs + vol_count]
            fine_face.e_g = e[N_surfs + vol_count]
            fine_face.r_g = r[N_surfs + vol_count]
            fine_face.g_g = g[N_surfs + vol_count]
            fine_face.q_g = Q[N_surfs + vol_count]

            for (wall_idx, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surf_count += 1
                    fine_face.T_w[wall_idx] = T[surf_count]
                    fine_face.j_w[wall_idx] = j[surf_count]
                    fine_face.g_a_w[wall_idx] = Abs[surf_count]
                    fine_face.e_w[wall_idx] = e[surf_count]
                    fine_face.r_w[wall_idx] = r[surf_count]
                    fine_face.g_w[wall_idx] = g[surf_count]
                    fine_face.q_w[wall_idx] = Q[surf_count]
                end
            end
        end
    end
    println("Done writing results back into the mesh.")

    return B, S_infty, A, R, C, D, M
end

# CLEAN Memory-optimized steady state - 2 matrices, no duplicates, single solve

# Minimal workspace with 2 matrices
mutable struct TwoMatrixWorkspace{P}
    matrix1::Matrix{P}    # B → M (system matrix)
    matrix2::Matrix{P}    # K → S_infty (preserved)
    
    # Working vectors
    work_vec1::Vector{P}
    work_vec2::Vector{P}
    work_vec3::Vector{P}
    work_vec4::Vector{P}
    work_vec5::Vector{P}
    
    # Physical property vectors
    Area::Vector{P}
    epsw::Vector{P}
    Tw::Vector{P}
    qw::Vector{P}
    bin_Qw_known::Vector{Int}
    Volume::Vector{P}
    Tg::Vector{P}
    qg::Vector{P}
    bin_Qg_known::Vector{Int}
    
    # Constructor
    function TwoMatrixWorkspace{P}(N_surfs::Int, N_vols::Int) where P
        n = N_surfs + N_vols
        new{P}(
            zeros(P, n, n), zeros(P, n, n),  # 2 matrices
            zeros(P, n), zeros(P, n), zeros(P, n), zeros(P, n), zeros(P, n),  # 5 vectors
            zeros(P, N_surfs), zeros(P, N_surfs), zeros(P, N_surfs), zeros(P, N_surfs), zeros(Int, N_surfs),
            zeros(P, N_vols), zeros(P, N_vols), zeros(P, N_vols), zeros(Int, N_vols)
        )
    end
end

# Populate workspace arrays from mesh
function populate_workspace!(ws::TwoMatrixWorkspace{P}, mesh::RayTracingMeshOptim) where P
    surf_count = 0
    vol_count = 0
    
    for i in 1:length(mesh.coarse_mesh)
        for (k, fine_face) in enumerate(mesh.fine_mesh[i])
            vol_count += 1
            ws.Volume[vol_count] = fine_face.volume
            ws.Tg[vol_count] = fine_face.T_in_g
            ws.qg[vol_count] = fine_face.q_in_g
            ws.bin_Qg_known[vol_count] = fine_face.T_in_g < 0.0 ? 1 : 0

            for (wall_idx, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surf_count += 1
                    ws.Area[surf_count] = fine_face.area[wall_idx]
                    ws.epsw[surf_count] = fine_face.epsilon[wall_idx]
                    ws.Tw[surf_count] = fine_face.T_in_w[wall_idx]
                    ws.qw[surf_count] = fine_face.q_in_w[wall_idx]
                    ws.bin_Qw_known[surf_count] = fine_face.T_in_w[wall_idx] < 0.0 ? 1 : 0
                end
            end
        end
    end
end

# Compute emissive powers
function compute_emissive_powers!(E_known::Vector{P}, Q_known::Vector{P}, ws::TwoMatrixWorkspace{P}, 
                                gas::GasProperties, N_surfs::Int, N_vols::Int) where P
    sigma = 5.670374419e-8
    
    # Build Q_vec_known
    Q_vec_known = vcat(ws.bin_Qw_known, ws.bin_Qg_known)
    
    # Surface emissive powers
    for i in 1:N_surfs
        if Q_vec_known[i] == 0
            E_known[i] = ws.epsw[i] * sigma * ws.Area[i] * ws.Tw[i]^4
            Q_known[i] = 0
        else
            E_known[i] = 0
            Q_known[i] = ws.qw[i] * ws.Area[i]
        end
    end
    
    # Volume emissive powers
    for i in 1:N_vols
        idx = N_surfs + i
        if Q_vec_known[idx] == 0
            E_known[idx] = 4 * gas.kappa * sigma * ws.Volume[i] * ws.Tg[i]^4
            Q_known[idx] = 0
        else
            E_known[idx] = 0
            Q_known[idx] = ws.qg[i] * ws.Volume[i]
        end
    end
end

# Write results back to mesh
function write_results_to_mesh!(mesh::RayTracingMeshOptim, T::Vector{P}, j::Vector{P}, 
                               Abs::Vector{P}, r::Vector{P}, N_surfs::Int) where P
    surf_count = 0
    vol_count = 0
    
    for (i, coarse_face) in enumerate(mesh.coarse_mesh)
        for (k, fine_face) in enumerate(mesh.fine_mesh[i])
            vol_count += 1
            vol_idx = N_surfs + vol_count
            
            # Volume properties
            e_vol = max(j[vol_idx] - r[vol_idx], 0.0)
            fine_face.T_g = T[vol_idx]
            fine_face.j_g = j[vol_idx]
            fine_face.g_a_g = Abs[vol_idx]
            fine_face.e_g = e_vol
            fine_face.r_g = r[vol_idx]
            fine_face.g_g = Abs[vol_idx] + r[vol_idx]
            fine_face.q_g = e_vol - Abs[vol_idx]

            # Surface properties
            for (wall_idx, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surf_count += 1
                    e_surf = max(j[surf_count] - r[surf_count], 0.0)
                    
                    fine_face.T_w[wall_idx] = T[surf_count]
                    fine_face.j_w[wall_idx] = j[surf_count]
                    fine_face.g_a_w[wall_idx] = Abs[surf_count]
                    fine_face.e_w[wall_idx] = e_surf
                    fine_face.r_w[wall_idx] = r[surf_count]
                    fine_face.g_w[wall_idx] = Abs[surf_count] + r[surf_count]
                    fine_face.q_w[wall_idx] = e_surf - Abs[surf_count]
                end
            end
        end
    end
    
    println("Results written: $vol_count volumes, $surf_count surfaces")
end

# Main memory-optimized steady state function - CLEAN VERSION
function steadyState_optimized!(mesh::RayTracingMeshOptim, F::Matrix{P}, gas::GasProperties) where {P<:Real}
    println("=== Memory-Optimized Steady State Solver ===")
    
    # Count surfaces and volumes
    N_surfs = sum([sum(fine_face.solidWalls) for i in 1:length(mesh.coarse_mesh) for fine_face in mesh.fine_mesh[i]])
    N_vols = length([fine_face for i in 1:length(mesh.coarse_mesh) for fine_face in mesh.fine_mesh[i]])
    n = N_surfs + N_vols
    
    println("Found $N_surfs surfaces and $N_vols volumes (total: $n elements).")
    
    # Check memory requirements
    matrix_gb = 2 * sizeof(P) * n^2 / 1e9
    total_gb = matrix_gb + sizeof(P) * n * 10 / 1e9
    println("Memory needed: ~$(total_gb) GB (2 matrices: $(matrix_gb) GB + vectors)")
    
    # Allocate workspace
    println("Allocating workspace...")
    ws = TwoMatrixWorkspace{P}(N_surfs, N_vols)
    
    # Use work vectors for main arrays
    E_known = ws.work_vec1
    Q_known = ws.work_vec2
    
    # Step 1: Populate workspace from mesh
    println("Populating workspace from mesh...")
    populate_workspace!(ws, mesh)
    
    # Step 2: Compute emissive powers
    println("Computing emissive powers...")
    compute_emissive_powers!(E_known, Q_known, ws, gas, N_surfs, N_vols)
    
    # Step 3: Compute B matrix in matrix1
    println("Computing B matrix...")
    B = ws.matrix1
    fill!(B, zero(P))
    
    if gas.omega > 1e-6 || sum(ws.epsw) < n
        @threads for i in 1:n
            @inbounds for j in 1:N_surfs
                B[i, j] = 1.0 - ws.epsw[j]
            end
            @inbounds for j in (N_surfs+1):n
                B[i, j] = gas.omega
            end
        end
    end
    
    # Step 4: Compute K matrix in matrix2
    println("Computing K matrix...")
    K = ws.matrix2
    @threads for i in 1:n
        @inbounds for j in 1:n
            if B[i, j] != 0.0 && F[i, j] != 0.0
                K[i, j] = F[i, j] * B[i, j]
            else
                K[i, j] = 0.0
            end
        end
    end
    
    # Step 5: Solve for S_infty (overwrite K in matrix2)
    println("Solving for S_infty...")
    if isapprox(gas.sigma_s, 0.0, atol=1e-10)
        copyto!(ws.matrix2, F)  # S_infty = F
    else
        # Convert K to (I - K) in matrix2
        @threads for i in 1:n
            @inbounds for j in 1:n
                if i == j
                    K[i, j] = 1.0 - K[i, j]
                else
                    K[i, j] = -K[i, j]
                end
            end
        end
        ws.matrix2 .= K \ F  # Solve and store in matrix2
    end
    
    S_infty = ws.matrix2  # matrix2 now contains S_infty and we keep it
    
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
    
    # Build system matrix M by computing A and R elements on-the-fly
    @threads for i in 1:n
        @inbounds for j in 1:n
            # Compute B values on-the-fly
            B_ji = if j <= N_surfs
                1.0 - ws.epsw[j]
            else
                gas.omega
            end
            
            B_ij = if i <= N_surfs
                1.0 - ws.epsw[i]
            else
                gas.omega
            end
            
            # A[i,j] = (1 - B[j,i]) * S_infty[i,j] * (1 - B[i,j])
            # R[i,j] = (1 - B[j,i]) * S_infty[i,j] * B[i,j]
            common_term = (1.0 - B_ji) * S_infty[i, j]
            A_ij = common_term * (1.0 - B_ij)
            R_ij = common_term * B_ij
            
            # Build system matrix based on boundary condition
            if Q_vec_known[i] == 1
                # M[i,j] = I[i,j] - A'[i,j] - R'[i,j] = I[i,j] - A[j,i] - R[j,i]
                if i == j
                    M[i, j] = 1.0 - A_ij - R_ij
                else
                    # Need A[j,i] and R[j,i] - compute with swapped indices
                    B_ij_swap = if j <= N_surfs
                        1.0 - ws.epsw[j]
                    else
                        gas.omega
                    end
                    B_ji_swap = if i <= N_surfs
                        1.0 - ws.epsw[i]
                    else
                        gas.omega
                    end
                    common_term_swap = (1.0 - B_ji_swap) * S_infty[j, i]
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
                    B_ij_swap = if j <= N_surfs
                        1.0 - ws.epsw[j]
                    else
                        gas.omega
                    end
                    B_ji_swap = if i <= N_surfs
                        1.0 - ws.epsw[i]
                    else
                        gas.omega
                    end
                    common_term_swap = (1.0 - B_ji_swap) * S_infty[j, i]
                    R_ji = common_term_swap * B_ij_swap
                    M[i, j] = -R_ji
                end
            end
        end
    end
    
    # Step 7: Solve linear system (SINGLE SOLVE)
    println("Solving linear system...")
    j = M \ b
    
    # Step 8: Compute Abs = A' * j and r = R' * j
    println("Computing absorbed and reflected energies...")
    Abs = ws.work_vec4
    r = ws.work_vec5
    
    fill!(Abs, zero(P))
    fill!(r, zero(P))
    
    # Compute matrix-vector products by recomputing A and R elements
    @threads for i in 1:n
        local_abs = zero(P)
        local_r = zero(P)
        
        @inbounds for k in 1:n
            # Compute A[k,i] and R[k,i]
            B_ik = if i <= N_surfs
                1.0 - ws.epsw[i]
            else
                gas.omega
            end
            B_ki = if k <= N_surfs
                1.0 - ws.epsw[k]
            else
                gas.omega
            end
            
            common_term = (1.0 - B_ik) * S_infty[k, i]
            A_ki = common_term * (1.0 - B_ki)
            R_ki = common_term * B_ki
            
            local_abs += A_ki * j[k]
            local_r += R_ki * j[k]
        end
        
        Abs[i] = local_abs
        r[i] = local_r
    end
    
    # Step 9: Compute temperatures
    println("Computing temperatures...")
    T = ws.work_vec1  # Reuse work_vec1 (E_known no longer needed)
    sigma = 5.670374419e-8
    
    @threads for i in 1:n
        if i <= N_surfs
            e_i = max(j[i] - r[i], 0.0)
            T[i] = (e_i / (ws.epsw[i] * sigma * ws.Area[i]))^0.25
        else
            vol_idx = i - N_surfs
            e_i = max(j[i] - r[i], 0.0)
            T[i] = (e_i / (4 * gas.kappa * ws.Volume[vol_idx] * sigma))^0.25
        end
        if isnan(T[i])
            T[i] = 0.0
        end
    end
    
    # Step 10: Write results to mesh
    println("Writing results to mesh...")
    write_results_to_mesh!(mesh, T, j, Abs, r, N_surfs)
    
    println("=== Steady State Solution Complete ===")
    println("Total memory used: ~$(total_gb) GB (2 matrices + vectors)")
    
    return nothing
end

function steadyState2D!(mesh::RayTracingMeshOptim, F::Matrix{P}, gas::GasProperties; return_matrices=false) where {P<:Real}
    if return_matrices
        return steadyState_simple!(mesh, F, gas)
    else
        return steadyState_optimized!(mesh, F, gas)    
    end
end