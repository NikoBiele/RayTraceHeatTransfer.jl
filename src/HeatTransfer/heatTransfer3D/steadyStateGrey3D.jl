# Simplified workspace for 3D surfaces (no volumes!)
mutable struct SurfaceOnlyWorkspace{P}
    matrix1::Matrix{P}    # B → M (system matrix)
    matrix2::Matrix{P}    # K → S_infty (preserved)
    
    # Working vectors
    work_vec1::Vector{P}
    work_vec2::Vector{P}
    work_vec3::Vector{P}
    work_vec4::Vector{P}
    work_vec5::Vector{P}
    
    # Physical property vectors (surfaces only!)
    Area::Vector{P}
    epsw::Vector{P}
    Tw::Vector{P}
    qw::Vector{P}
    bin_Qw_known::Vector{Int}
    
    # Constructor
    function SurfaceOnlyWorkspace{P}(N_surfs::Int) where P
        new{P}(
            zeros(P, N_surfs, N_surfs), zeros(P, N_surfs, N_surfs),  # 2 matrices
            zeros(P, N_surfs), zeros(P, N_surfs), zeros(P, N_surfs), 
            zeros(P, N_surfs), zeros(P, N_surfs),  # 5 vectors
            zeros(P, N_surfs), zeros(P, N_surfs), zeros(P, N_surfs), 
            zeros(P, N_surfs), zeros(Int, N_surfs)
        )
    end
end

# Populate workspace from 3D domain (surfaces only!)
function populate_workspace_3D!(ws::SurfaceOnlyWorkspace{P}, domain::Domain3D_faces, 
                                spectral_bin::Int) where P
    surf_count = 0
    
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            ws.Area[surf_count] = subface.area
            
            # Handle spectral or grey epsilon
            if isa(subface.epsilon, Vector)
                ws.epsw[surf_count] = subface.epsilon[spectral_bin]
            else
                ws.epsw[surf_count] = subface.epsilon
            end
            
            ws.Tw[surf_count] = subface.T_in_w
            ws.qw[surf_count] = subface.q_in_w
            ws.bin_Qw_known[surf_count] = subface.T_in_w < 0.0 ? 1 : 0
        end
    end
end

# Compute emissive powers (surfaces only!)
function compute_emissive_powers_3D!(E_known::Vector{P}, Q_known::Vector{P}, 
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
function compute_temperatures_3D!(T::Vector{P}, j::Vector{P}, r::Vector{P}, 
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

# Helper to get B value for surfaces (always just reflectivity)
@inline function get_B_value_3D(ws::SurfaceOnlyWorkspace{P}, idx::Int) where {P}
    return 1.0 - ws.epsw[idx]
end

# Main grey solver for 3D surfaces - ALWAYS solves
function steadyStateGrey3D!(domain::Domain3D_faces, F::Matrix{P}; 
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
    populate_workspace_3D!(ws, domain, spectral_bin)
    
    # Step 2: Compute emissive powers
    println("Computing emissive powers...")
    compute_emissive_powers_3D!(E_known, Q_known, ws, N_surfs)
    
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
    
    S_infty = ws.matrix2
    
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
            B_ji = get_B_value_3D(ws, j)
            B_ij = get_B_value_3D(ws, i)
            
            common_term = (1.0 - B_ji) * S_infty[i, j]
            A_ij = common_term * (1.0 - B_ij)
            R_ij = common_term * B_ij
            
            if Q_vec_known[i] == 1  # Known flux
                # Need A[j,i] and R[j,i]
                B_ij_swap = get_B_value_3D(ws, j)
                B_ji_swap = get_B_value_3D(ws, i)
                common_term_swap = (1.0 - B_ji_swap) * S_infty[j, i]
                A_ji = common_term_swap * (1.0 - B_ij_swap)
                R_ji = common_term_swap * B_ij_swap
                
                M[i, j] = (i == j) ? (1.0 - A_ij - R_ij) : (-A_ji - R_ji)
            else  # Known temperature
                # Need R[j,i]
                B_ij_swap = get_B_value_3D(ws, j)
                B_ji_swap = get_B_value_3D(ws, i)
                common_term_swap = (1.0 - B_ji_swap) * S_infty[j, i]
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
        
        B_i = get_B_value_3D(ws, i)
        
        for k in 1:n
            B_k = get_B_value_3D(ws, k)
            
            common_term = (1.0 - B_i) * S_infty[k, i]
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
    compute_temperatures_3D!(T, j, r, ws, N_surfs)
    
    # Step 10: Write results to domain
    println("Writing results to domain...")
    write_results_to_domain_3D!(domain, T, j, Abs, r; spectral_bin=spectral_bin)

    # Step 11: Compute energy conservation error
    println("Computing energy conservation error...")
    domain.energy_error = sum(j - r - Abs)
    
    println("=== 3D Grey Solution Complete ===")
    
    return nothing
end

function build_system_matrices3D!(domain::Domain3D_faces, F::Matrix{P}; 
                                   spectral_bin::Int=1) where {P<:Real}
    
    # Count surfaces (no volumes in 3D!)
    N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    n = N_surfs
    
    # Allocate output matrices
    B_save = zeros(P, n, n)
    S_infty_save = zeros(P, n, n)
    A_save = zeros(P, n, n)
    R_save = zeros(P, n, n)
    C_save = zeros(P, n, n)
    D_save = zeros(P, n, n)
    M_save = zeros(P, n, n)
    
    # Allocate workspace
    ws = SurfaceOnlyWorkspace{P}(N_surfs)
    
    # Use work vectors for emissive powers (needed for M matrix)
    E_known = ws.work_vec1
    Q_known = ws.work_vec2
    
    # Step 1: Populate workspace from domain
    populate_workspace_3D!(ws, domain, spectral_bin)
    
    # Step 2: Compute emissive powers (needed for M matrix boundary conditions)
    compute_emissive_powers_3D!(E_known, Q_known, ws, N_surfs)
    
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
    B_save .= B
    
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
    
    S_infty = ws.matrix2  # matrix2 now contains S_infty
    S_infty_save .= S_infty
    
    # Step 6: Compute A (absorption) and R (reflection) matrices
    for i in 1:n
        for j in 1:n
            # Get B values using helper function
            B_ji = get_B_value_3D(ws, j)
            B_ij = get_B_value_3D(ws, i)
            
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
    Q_vec_known = ws.bin_Qw_known
    
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