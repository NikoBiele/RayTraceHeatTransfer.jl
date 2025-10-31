# Fixed radiative transfer algorithm
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
    kappa_g::Vector{P}      # NEW: absorption coefficient per volume
    omega_g::Vector{P}      # NEW: scattering albedo per volume
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
            zeros(P, N_vols), zeros(P, N_vols), zeros(P, N_vols), zeros(P, N_vols), zeros(P, N_vols), zeros(Int, N_vols)
        )
    end
end

# Populate workspace arrays from mesh - NOW with kappa_g and omega_g precomputation
function populate_workspace!(ws::TwoMatrixWorkspace{P}, mesh::RayTracingDomain2D, spectral_bin::Int) where P
    surf_count = 0
    vol_count = 0
    
    for i in 1:length(mesh.coarse_mesh)
        for (k, fine_face) in enumerate(mesh.fine_mesh[i])
            vol_count += 1
            ws.Volume[vol_count] = fine_face.volume
            
            # NEW: Store kappa and precompute omega
            ws.kappa_g[vol_count] = fine_face.kappa_g[spectral_bin]
            local_kappa = fine_face.kappa_g[spectral_bin]
            local_sigma_s = fine_face.sigma_s_g[spectral_bin]
            local_beta = local_kappa + local_sigma_s
            ws.omega_g[vol_count] = local_beta > 0.0 ? local_sigma_s / local_beta : 0.0
            
            ws.Tg[vol_count] = fine_face.T_in_g
            ws.qg[vol_count] = fine_face.q_in_g
            ws.bin_Qg_known[vol_count] = fine_face.T_in_g < 0.0 ? 1 : 0

            for (wall_idx, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surf_count += 1
                    ws.Area[surf_count] = fine_face.area[wall_idx]
                    ws.epsw[surf_count] = fine_face.epsilon[wall_idx][spectral_bin]
                    ws.Tw[surf_count] = fine_face.T_in_w[wall_idx]
                    ws.qw[surf_count] = fine_face.q_in_w[wall_idx]
                    ws.bin_Qw_known[surf_count] = fine_face.T_in_w[wall_idx] < 0.0 ? 1 : 0
                end
            end
        end
    end
end

# Updated compute_emissive_powers - use precomputed kappa
function compute_emissive_powers_variable!(E_known::Vector{P}, Q_known::Vector{P}, ws::TwoMatrixWorkspace{P}, 
                                          mesh::RayTracingDomain2D,
                                          N_surfs::Int, N_vols::Int) where {P}
    
    # Build Q_vec_known
    Q_vec_known = vcat(ws.bin_Qw_known, ws.bin_Qg_known)
    
    # Surface emissive powers
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
    for vol_count in 1:N_vols
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
function compute_temperatures_variable!(T::Vector{P}, j::Vector{P}, r::Vector{P}, ws::TwoMatrixWorkspace{P},
                                       mesh::RayTracingDomain2D,
                                       N_surfs::Int, N_vols::Int, spectral_bin::Int) where {P}
    
    # Surface temperatures
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
    for vol_count in 1:N_vols
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

# NEW: Helper to get B value using precomputed omega
@inline function get_B_value(ws::TwoMatrixWorkspace{P}, idx::Int, N_surfs::Int) where {P}
    if idx <= N_surfs
        return 1.0 - ws.epsw[idx]
    else
        return ws.omega_g[idx - N_surfs]
    end
end

# Main solver - UNCHANGED except using get_B_value
function steadyStateGrey2D!(mesh::RayTracingDomain2D, F::Matrix{P}; spectral_bin::N=1) where {P<:Real, N<:Integer}
    println("=== Variable Extinction Memory-Optimized Steady State Solver ===")
    
    # Count surfaces and volumes
    N_surfs = sum([sum(fine_face.solidWalls) for i in 1:length(mesh.coarse_mesh) for fine_face in mesh.fine_mesh[i]])
    N_vols = length([fine_face for i in 1:length(mesh.coarse_mesh) for fine_face in mesh.fine_mesh[i]])
    n = N_surfs + N_vols
    
    println("Found $N_surfs surfaces and $N_vols volumes (total: $n elements).")
    
    # Allocate workspace
    println("Allocating workspace...")
    ws = TwoMatrixWorkspace{P}(N_surfs, N_vols)
    
    # Use work vectors for main arrays
    E_known = ws.work_vec1
    Q_known = ws.work_vec2
    
    # Step 1: Populate workspace from mesh
    println("Populating workspace from mesh...")
    populate_workspace!(ws, mesh, spectral_bin)
    
    # Step 2: Compute emissive powers with variable extinction
    println("Computing emissive powers with variable extinction...")
    compute_emissive_powers_variable!(E_known, Q_known, ws, mesh, N_surfs, N_vols)
    
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
    
    S_infty = ws.matrix2  # matrix2 now contains S_infty
    
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
    
    # Build system matrix M - NOW using get_B_value instead of get_local_omega
    for i in 1:n
        for j in 1:n
            # Get B values using helper (no nested loops!)
            B_ji = get_B_value(ws, j, N_surfs)
            B_ij = get_B_value(ws, i, N_surfs)
            
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
                    B_ij_swap = get_B_value(ws, j, N_surfs)
                    B_ji_swap = get_B_value(ws, i, N_surfs)
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
                    B_ij_swap = get_B_value(ws, j, N_surfs)
                    B_ji_swap = get_B_value(ws, i, N_surfs)
                    common_term_swap = (1.0 - B_ji_swap) * S_infty[j, i]
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
    
    # Compute matrix-vector products - NOW using get_B_value
    for i in 1:n
        local_abs = zero(P)
        local_r = zero(P)
        
        for k in 1:n
            # Get B values using helper
            B_ik = get_B_value(ws, i, N_surfs)
            B_ki = get_B_value(ws, k, N_surfs)
            
            common_term = (1.0 - B_ik) * S_infty[k, i]
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
    compute_temperatures_variable!(T, j, r, ws, mesh, N_surfs, N_vols, 1)
    
    # Step 10: Write results to mesh
    println("Writing results to mesh...")
    write_results_to_mesh!(mesh, T, j, Abs, r, N_surfs; spectral_bin=spectral_bin)

    # Step 11: Compute energy conservation error
    println("Computing energy conservation error...")
    mesh.energy_error = sum(j - r - Abs)
    
    println("=== Variable Extinction Steady State Solution Complete ===")
    
    return nothing
end

function build_system_matrices2D!(mesh::RayTracingDomain2D, F::Matrix{P}; 
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
    populate_workspace!(ws, mesh, spectral_bin)
    
    # Step 2: Compute emissive powers (needed for M matrix boundary conditions)
    compute_emissive_powers_variable!(E_known, Q_known, ws, mesh, N_surfs, N_vols)
    
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
            B_ji = get_B_value(ws, j, N_surfs)
            B_ij = get_B_value(ws, i, N_surfs)
            
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
    
    return B_save, S_infty_save, A_save, R_save, C_save, D_save, M_save
end