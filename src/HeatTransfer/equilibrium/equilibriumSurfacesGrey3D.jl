
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
function equilibriumSurfacesGrey3D!(domain::ViewFactorDomain3D, F::AbstractMatrix; 
                                   spectral_bin::Int=1, verbose::Bool=true) #where {P<:Real}
    verbose && println("=== 3D Surface-Only Grey Solver ===")
    
    P = Float64

    # Count surfaces (no volumes!)
    N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    n = N_surfs
    
    verbose && println("Found $N_surfs surfaces")
    
    # Allocate workspace
    ws = SurfaceOnlyWorkspace{P}(N_surfs)
    
    # Use work vectors
    E_known = ws.work_vec1
    Q_known = ws.work_vec2
    
    # Populate workspace
    verbose && println("Populating workspace...")
    populateWorkspace!(ws, domain, spectral_bin)
    
    # Compute emissive powers
    verbose && println("Computing emissive powers...")
    computeEmissivePowers!(E_known, Q_known, ws, N_surfs)
    
    # Compute B matrix (surface reflectivity only!)
    verbose && println("Computing B matrix...")
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
    
    # Build system matrix M
    verbose && println("Assembling linear system...")
    M = ws.matrix1  # Reuse matrix1
    h = ws.work_vec3
    
    Q_vec_known = ws.bin_Qw_known
    for i in 1:n
        h[i] = Q_vec_known[i] == 1 ? Q_known[i] : E_known[i]
    end
    
    # Build system matrix M  ==  I − Diagonal(coeff)·Fᵀ
    coeff = ifelse.(Q_vec_known .== 1, 1.0, 1.0 .- ws.epsw)
    M .= I - Diagonal(coeff) * permutedims(F)
    
    # Step 7: Solve
    verbose && println("Solving linear system...")
    j = M \ h
    
    # Step 8: Compute absorbed and reflected
    verbose && println("Computing absorbed and reflected energies...")
    Abs = ws.work_vec4
    r = ws.work_vec5
    fill!(Abs, zero(P))
    fill!(r, zero(P))
    
    for i in 1:n
        B_i = getValueB(ws, i)
        g_i = zero(P)
        for k in 1:n
            g_i += F[k, i] * j[k]
        end
        Abs[i] = (1.0 - B_i) * g_i
        r[i]   = B_i * g_i
    end
    
    # Step 9: Compute temperatures
    verbose && println("Computing temperatures...")
    T = ws.work_vec1
    computeTemperatures!(T, j, r, ws, N_surfs)
    
    # Step 10: Write results to domain
    verbose && println("Writing results to domain...")
    writeResultsToDomain!(domain, j, Abs, r; T=T, spectral_bin=spectral_bin)

    # Step 11: Compute energy conservation error
    verbose && println("Computing energy conservation error...")
    domain.energy_error = sum(j - r - Abs)
    
    verbose && println("=== 3D Grey Solution Complete ===")
    
    return nothing
end

