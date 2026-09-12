# Hooks used by the spectral solvers — thin dispatch onto the spectral model.

"""
    validateSpectralSetup(domain)

Check that a spectral domain has a `spectral_model` consistent with its
number of bins. Called at the entry of every spectral solver.
"""
function validateSpectralSetup(domain)
    m = domain.spectral_model
    if m === nothing
        error("""
        Spectral solve requires a spectral model to be set, e.g.

            mesh.spectral_model = PlanckBands(10 .^ range(log10(1e-7), log10(1e-3), length=51))
            mesh.spectral_model = adaptiveSpectralBins(λ, κ; tol=1e-3, L_range=(h, L), T_range=(T_lo, T_hi))
            mesh.spectral_model = ConstantWeights(a)
        """)
    end
    validate(m, domain.n_spectral_bins)
    return nothing
end

"""
    _fill_bin_fractions!(emitFrac, i, T, domain)

Row `i` of the per-bin emission fractions at temperature `T`.
"""
function _fill_bin_fractions!(emitFrac::AbstractMatrix, i::Int, T::Real, domain)
    fill_bin_fractions!(view(emitFrac, i, :), domain.spectral_model, T)
    return nothing
end

"""
    getBinsEmissionFractionDerivatives(domain, temperatures, N) -> Matrix (N × K)

Φ[i,k] = ∂f_k/∂T at the temperature of element i, for the first `N` elements.
"""
function getBinsEmissionFractionDerivatives(domain, temperatures::AbstractVector{G}, N::Int) where {G}
    Phi = zeros(G, N, domain.n_spectral_bins)
    m = domain.spectral_model
    for i in 1:N
        bin_fraction_derivatives!(view(Phi, i, :), m, temperatures[i])
    end
    return Phi
end