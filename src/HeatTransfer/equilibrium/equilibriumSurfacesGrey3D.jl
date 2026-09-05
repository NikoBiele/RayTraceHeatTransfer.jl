
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

"""
    SurfaceSystemOp(F, coeff, n)

Matrix-free form of `M == I − Diagonal(coeff)·Fᵀ`.

`M·x = x − coeff .* (Fᵀx)`, and `Fᵀx` is the cache-friendly direction for a
CSC matrix — each output entry is a dot product of one column with `x`, so no
transpose is materialised. Forming `M` explicitly would allocate three
full-size sparse temporaries.
"""
struct SurfaceSystemOp{T, MT}
    F::MT
    coeff::Vector{T}
    n::Int
    tmp::Vector{T}
end

SurfaceSystemOp(F, coeff::Vector{T}, n::Int) where {T} =
    SurfaceSystemOp{T, typeof(F)}(F, coeff, n, Vector{T}(undef, n))

Base.size(A::SurfaceSystemOp)                = (A.n, A.n)
Base.size(A::SurfaceSystemOp, ::Integer)     = A.n
Base.eltype(::SurfaceSystemOp{T}) where {T}  = T

function LinearAlgebra.mul!(y::AbstractVector, A::SurfaceSystemOp, x::AbstractVector)
    mul!(y, transpose(A.F), x)
    @inbounds @. y = x - A.coeff * y
    return y
end

function LinearAlgebra.mul!(y::AbstractVector, A::SurfaceSystemOp, x::AbstractVector,
                            α::Number, β::Number)
    mul!(A.tmp, transpose(A.F), x)
    @inbounds @. A.tmp = x - A.coeff * A.tmp
    @inbounds @. y = α * A.tmp + β * y
    return y
end

# Main grey solver for 3D surfaces - ALWAYS solves
function equilibriumSurfacesGrey3D!(domain::D, F::AbstractMatrix; 
                                   spectral_bin::Int=1, verbose::Bool=true) where {D<:SurfaceDomain3D}
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

    # Step 3: Compute B matrix with variable extinction
    verbose && println("Computing B matrix with variable extinction...")
    b = zeros(n)
    
    # Check if we need scattering calculations
    has_reflection = sum(ws.epsw) < n
    
    if has_reflection
        # Surface reflectivity terms
        for j in 1:N_surfs
            b[j] = 1.0 - ws.epsw[j]
        end
    end
    
    # Build system matrix M
    verbose && println("Assembling linear system...")
    h = ws.work_vec3
    
    Q_vec_known = ws.bin_Qw_known
    for i in 1:n
        h[i] = Q_vec_known[i] == 1 ? Q_known[i] : E_known[i]
    end
    
    # Build system matrix M  ==  I − Diagonal(coeff)·Fᵀ
    coeff = ifelse.(Q_vec_known .== 1, 1.0, 1.0 .- ws.epsw)

    # Step 7: Solve linear system
    verbose && println("Solving linear system...")
    if F isa SparseMatrixCSC
        op  = SurfaceSystemOp(F, coeff, n)
        wkr = GmresWorkspace(op, h; memory = 50)  # GMRES(50): hard cap on the Krylov subspace
        gmres!(wkr, op, h; restart = true, rtol = 1e-12)
        j = wkr.x          # solution; wkr.stats for residual history / iters
    else
        M = ws.matrix1
        M .= I - Diagonal(coeff) * F'
        j = M \ h
    end

    # Step 6: Compute Abs = A' * j and r = R' * j
    verbose && println("Computing absorbed and reflected energies...")
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

