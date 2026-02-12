function chooseMatrixType(N_elements::P, n_spectral_bins::P, F::Union{Matrix{G}, Vector{Matrix{G}}}) where {P<:Integer, G}
    total_size = N_elements * n_spectral_bins
    
    # Dense memory estimate (GB)
    dense_gb = 2 * total_size^2 * sizeof(G) / 1e9
    
    # Choose sparse or dense
    if G == BigFloat && dense_gb < 4.0
        println("Using high precision: Dense memory estimate = $(round(dense_gb, digits=2)) GB")
        return :dense  # High precision if fits in memory
    elseif dense_gb > 1.0  # 1 GB threshold
        return :sparse
    elseif total_size > 5000
        return :sparse
    else
        return :dense
    end
end