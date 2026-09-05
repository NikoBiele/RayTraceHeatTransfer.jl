"""
    _blackbody_F(λ, T)

Cumulative blackbody fraction F(0→λT) at a single wavelength λ [m] and
temperature T [K]. Thin wrapper so channel code can evaluate at arbitrary
wavelengths; reuses the existing series implementation.
"""
_blackbody_F(λ::Real, T::Real) = emitFracBlackBodySpectrum((λ,), T, 1)

"""
    _fill_emission_fractions!(emitFrac, i, T, limits, n_bins)

Bin-mode row fill: bin k spans [limits[k], limits[k+1]]; first bin extended to
0, last bin extended to ∞ (fractions sum to 1).
"""
function _fill_emission_fractions!(emitFrac::AbstractMatrix, i::Int, T::Real,
                                   limits::AbstractVector, n_bins::Int)
    F_prev = 0.0                                   # first bin extends down to λ=0
    for k in 1:n_bins
        if k == n_bins
            emitFrac[i, k] = 1.0 - F_prev          # last bin extends up to λ=∞
        else
            F_curr = emitFracBlackBodySpectrum(limits, T, k + 1)  # upper edge of bin k
            emitFrac[i, k] = F_curr - F_prev
            F_prev = F_curr
        end
    end
    return nothing
end

function getBinsEmissionFractions(rtm::RayTracingDomain2D, temperatures::Vector{G}) where {G}

    num_surfaces = length(rtm.surface_mapping)
    n_elem = num_surfaces + length(rtm.volume_mapping)
    emitFrac = zeros(G, n_elem, rtm.n_spectral_bins)

    for i in 1:n_elem
        _fill_emission_fractions!(emitFrac, i, temperatures[i],
                                      rtm.wavelength_band_limits, rtm.n_spectral_bins)
    end

    return emitFrac
end

function getBinsEmissionFractions(domain::SurfaceDomain3D{G,P}, temperatures::Vector{G}) where {G,P}

    num_surfaces = sum([length(superface.subFaces) for superface in domain.facesMesh])
    emitFrac = zeros(G, num_surfaces, domain.n_spectral_bins)

    for i in 1:num_surfaces
        _fill_emission_fractions!(emitFrac, i, temperatures[i],
                                    domain.wavelength_band_limits, domain.n_spectral_bins)
    end

    return emitFrac
end

"""
    getWeightedEmissionFractions(rtm, temperatures) -> Matrix (n_elem × K)

Per-band emission fractions weighted by the band property:
surfaces  f_k = ε_k·w_k(T) / Σ_i ε_i·w_i(T)
volumes   f_k = κ_k·w_k(T) / Σ_i κ_i·w_i(T)
Rows sum to 1. Reduces to plain Planck fractions when the property is
bin-uniform (then the property cancels).
"""
function getWeightedEmissionFractions(rtm::RayTracingDomain2D, temperatures)
    w = getBinsEmissionFractions(rtm, temperatures)   # plain Planck fractions
    K = rtm.n_spectral_bins
    num_s = length(rtm.surface_mapping)
    out = similar(w)
    for ((ci, fi, wi), si) in rtm.surface_mapping
        eps_k = rtm.fine_mesh[ci][fi].epsilon[wi]      # Vector over bins
        @views out[si, :] .= eps_k .* w[si, :]
        s = sum(@view out[si, :])
        s > 0 && (@views out[si, :] ./= s)
    end
    for ((ci, fi), vi) in rtm.volume_mapping
        kap_k = rtm.fine_mesh[ci][fi].kappa_g
        row = num_s + vi
        @views out[row, :] .= kap_k .* w[row, :]
        s = sum(@view out[row, :])
        s > 0 && (@views out[row, :] ./= s)
    end
    return out
end

function getWeightedEmissionFractions(domain::SurfaceDomain3D, temperatures)
    w = getBinsEmissionFractions(domain, temperatures)
    out = similar(w)
    surf_count = 0
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            eps_k = isa(subface.epsilon, Vector) ? subface.epsilon : fill(subface.epsilon, domain.n_spectral_bins)
            @views out[surf_count, :] .= eps_k .* w[surf_count, :]
            s = sum(@view out[surf_count, :])
            s > 0 && (@views out[surf_count, :] ./= s)
        end
    end
    return out
end