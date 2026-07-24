# Fixed radiative transfer algorithm

# Updated compute_emissive_powers - use precomputed kappa
function computeEmissivePowersVariable!(mesh::RayTracingDomain2D, ws::OneMatrixWorkspace{P},
                                        E_known::Vector{P}, Q_known::Vector{P}) where {P}
    
    # Build Q_vec_known
    if mesh.surfaces_only
        Q_vec_known = vcat(ws.bin_Qw_known)
    else
        Q_vec_known = vcat(ws.bin_Qw_known, ws.bin_Qg_known)
    end

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
    if !mesh.surfaces_only
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
end

# Updated compute temperatures - use precomputed kappa
function computeTemperaturesVariable!(mesh::RayTracingDomain2D, ws::OneMatrixWorkspace{P},
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
    if !mesh.surfaces_only
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
end

# Main solver - UNCHANGED except using getValueB
function equilibriumGrey2D!(mesh::RayTracingDomain2D, F::AbstractMatrix; spectral_bin::N=1) where {N<:Integer} # P<:Real, 
    println("=== Variable Extinction Memory-Optimized Steady State Solver ===")
    
    # Count surfaces and volumes
    N_surfs = length(mesh.surface_mapping)
    N_vols = mesh.surfaces_only ? 0 : length(mesh.volume_mapping)
    n = N_surfs + N_vols
        
    # Allocate workspace
    println("Allocating workspace...")
    ws = OneMatrixWorkspace{Float64}(N_surfs, N_vols)
    
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
    b = zeros(n)
    
    # Check if we need scattering calculations
    has_scattering = any(ws.omega_g .> 1e-6)
    has_reflection = sum(ws.epsw) < n
    
    if has_scattering || has_reflection
        # Surface reflectivity terms
        for j in 1:N_surfs
            b[j] = 1.0 - ws.epsw[j]
        end
        
        # Volume scattering terms - use precomputed omega
        if !mesh.surfaces_only
            for vol_count in 1:N_vols
                vol_idx = N_surfs + vol_count
                b[vol_idx] = ws.omega_g[vol_count]
            end
        end
    end
    
    # Step 4: Build system matrix M directly in matrix1 (overwrite B)
    println("Assembling linear system...")
    M = ws.matrix1
    fill!(M, zero(Float64))
    h = zeros(n)
    
    # Build Q_vec_known and RHS
    if mesh.surfaces_only
        Q_vec_known = vcat(ws.bin_Qw_known)
    else
        Q_vec_known = vcat(ws.bin_Qw_known, ws.bin_Qg_known)
    end
    for i in 1:n
        if Q_vec_known[i] == 1
            h[i] = Q_known[i]
        else
            h[i] = E_known[i]
        end
    end
    
    # Build system matrix M  ==  I − Diagonal(coeff)·Fᵀ
    coeff = ifelse.(Q_vec_known .== 1, 1.0, b)
    M = I - Diagonal(coeff) * permutedims(F)
    
    # Step 5: Solve linear system
    println("Solving linear system...")
    if F isa SparseMatrixCSC
        wkr = GmresWorkspace(M, h; memory = 50)   # GMRES(50): hard cap on the Krylov subspace
        gmres!(wkr, M, h; restart = true, rtol = 1e-12)
        j = wkr.x          # solution; ws.stats for residual history / iters
    else
        j = M \ h
    end
        
    # Step 6: Compute Abs = A' * j and r = R' * j
    println("Computing absorbed and reflected energies...")
    Abs = ws.work_vec4
    r = ws.work_vec5
    fill!(Abs, zero(Float64))
    fill!(r, zero(Float64))
    
    if F isa SparseMatrixCSC
        # Column sweep over F (CSC-native): column i of F holds F[k,i] for all senders k.
        # g_i = Σ_k F[k,i]·j[k] is the total incident power on element i.
        # Receiver i then splits its OWN incident power:
        #   r[i]   = b[i]     · g_i     (reflected/scattered by receiver i)
        #   Abs[i] = (1−b[i]) · g_i     (absorbed by receiver i)
        rows = rowvals(F); vals = nonzeros(F)
        @inbounds for i in 1:n                 # column i == receiving element
            g_i = 0.0
            for idx in nzrange(F, i)
                k = rows[idx]                  # sending element
                g_i += vals[idx] * j[k]        # F[k,i]·j[k]
            end
            r[i]   = b[i] * g_i
            Abs[i] = (1.0 - b[i]) * g_i
        end
    else
        # Dense: same receiver-indexed split
        for i in 1:n
            g_i = zero(Float64)
            for k in 1:n
                g_i += F[k,i] * j[k]
            end
            r[i]   = b[i] * g_i
            Abs[i] = (1.0 - b[i]) * g_i
        end
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