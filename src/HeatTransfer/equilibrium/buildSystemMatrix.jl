function buildSystemMatrixVolume(mesh, F::AbstractMatrix; 
                                spectral_bin::Int=1) # where {P<:Real}
    
    # Count surfaces and volumes
    N_surfs = sum([sum(fine_face.solidWalls) for i in 1:length(mesh.coarse_mesh) 
                   for fine_face in mesh.fine_mesh[i]])
    N_vols = mesh.surfaces_only ? 0 : length([fine_face for i in 1:length(mesh.coarse_mesh) 
                     for fine_face in mesh.fine_mesh[i]])
    
    # Allocate workspace
    ws = OneMatrixWorkspace{Float64}(N_surfs, N_vols)
    
    # Use work vectors for emissive powers (needed for M matrix)
    E_known = ws.work_vec1
    Q_known = ws.work_vec2
    
    # Step 1: Populate workspace from mesh
    populateWorkspace!(ws, mesh, spectral_bin)
    
    # Step 2: Compute emissive powers (needed for M matrix boundary conditions)
    computeEmissivePowersVariable!(mesh, ws, E_known, Q_known)
    
    # Step 3: Compute B matrix (reflectivity/scattering albedo)
    b = get_b(mesh)

    # Build M matrix (system matrix for solving)
    M = ws.matrix1
    fill!(M, zero(Float64))
    Q_vec_known = vcat(ws.bin_Qw_known, ws.bin_Qg_known)
    
    # Build system matrix M  ==  I − Diagonal(coeff)·Fᵀ
    coeff = ifelse.(Q_vec_known .== 1, 1.0, b[:,spectral_bin])
    M = I - Diagonal(coeff) * F'
    
    return M
end

function buildSystemMatrixSurface(domain, F::AbstractMatrix; 
                                   spectral_bin::Int=1) #where {P<:Real}
    P = Float64

    # Count surfaces (no volumes in 3D!)
    if typeof(domain) <: ViewFactorDomain3D
        N_surfs = sum([length(superface.subFaces) for superface in domain.facesMesh])
    else
        N_surfs = length(domain.surface_mapping)
    end
    n = N_surfs
    
    # Allocate workspace
    ws = SurfaceOnlyWorkspace{Float64}(N_surfs)
    
    # Use work vectors for emissive powers (needed for M matrix)
    E_known = ws.work_vec1
    Q_known = ws.work_vec2
    
    # Step 1: Populate workspace from domain
    populateWorkspace!(ws, domain, spectral_bin)
    
    # Step 2: Compute emissive powers (needed for M matrix boundary conditions)
    computeEmissivePowers!(E_known, Q_known, ws, N_surfs)
    
    # Step 3: Compute B matrix (surface reflectivity only!)
    b = get_b(domain)
    
    # Build M matrix (system matrix for solving)
    Q_vec_known = ws.bin_Qw_known
    
    # Build system matrix M  ==  I − Diagonal(coeff)·Fᵀ
    coeff = ifelse.(Q_vec_known .== 1, 1.0, b[1:n,spectral_bin])
    M = I - Diagonal(coeff) * F'
    
    return M
end

function buildSystemMatrix(domain, F::AbstractMatrix;
                                   spectral_bin::Int=1) # where {P<:Real}
    if domain.surfaces_only
        return buildSystemMatrixSurface(domain, F; spectral_bin)
    else
        return buildSystemMatrixVolume(domain, F; spectral_bin)
    end
end