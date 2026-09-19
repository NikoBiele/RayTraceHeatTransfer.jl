# Reference solution for a grey plane-parallel slab in radiative equilibrium with
# anisotropic (Henyey–Greenstein) scattering, between two grey walls that reflect
# diffusely and/or specularly. Deterministic: 1D discrete ordinates, one linear solve.
#
#   * ordinates: double Gauss (Gauss–Legendre on each half-range), so the half-range
#     wall fluxes are integrated exactly;
#   * phase function: azimuthally averaged HG from its Legendre series
#     p(μ, μ′) = Σ_l (2l+1) gˡ P_l(μ) P_l(μ′), truncated at l = n_half − 1, the highest
#     order for which the discrete kernel is exactly normalised and symmetric;
#   * space: Nx cells with a cell-constant source and exact attenuation along every
#     ordinate — the same assumption as uniform re-emission within an element, so at
#     equal Nx the comparison with the package isolates the angular treatment;
#   * radiative equilibrium: absorption + isotropic re-emission is the isotropic part
#     of a conservative kernel, p_eff = (1 − ω) + ω p_HG.
#
# With g = 0 (or ω = 0) this is the grey slab of Heaslet & Warming (1965).
# Used by test/test_anisotropic_slab.jl and the README example; never loaded by the package.

using LinearAlgebra

const ASLAB_σ = 5.670374419e-8

# Gauss–Legendre nodes and weights on [−1, 1] (Golub–Welsch)
function aslab_gauss_legendre(n::Int)
    β = [k / sqrt(4.0 * k^2 - 1.0) for k in 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), β))
    return E.values, 2.0 .* E.vectors[1, :] .^ 2
end

# double-Gauss ordinates: entries 1:n_half have μ > 0, entries n_half+1:2n_half are their mirrors −μ
function aslab_ordinates(n_half::Int)
    x, w = aslab_gauss_legendre(n_half)
    μp = (x .+ 1.0) ./ 2.0
    wp = w ./ 2.0
    return vcat(μp, -μp), vcat(wp, wp)
end

# scattering matrix K[m, m′] = ½ w_m′ p_eff(μ_m, μ_m′): source S_m = Σ_m′ K[m, m′] Ī_m′
function aslab_kernel(μ::Vector{Float64}, w::Vector{Float64}, ω::Float64, g::Float64, Lmax::Int)
    n = length(μ)
    P = zeros(Lmax + 1, n)                           # P[l+1, m] = P_l(μ_m)
    P[1, :] .= 1.0
    Lmax >= 1 && (P[2, :] .= μ)
    for l in 1:Lmax-1
        @. P[l + 2, :] = ((2l + 1) * μ * P[l + 1, :] - l * P[l, :]) / (l + 1)
    end
    p = fill(1.0 - ω, n, n)
    for l in 0:Lmax
        p .+= ω * (2l + 1) * g^l .* (P[l + 1, :] * P[l + 1, :]')
    end
    return 0.5 .* p .* w'
end

"""
    anisotropic_slab_equilibrium(κ, σ_s, g, L, T1, T2; ε1 = 1, ε2 = 1, specularity = 0, Nx = 32, n_half = 32)

Temperatures of the `Nx` cells (from the wall at `T1` to the wall at `T2`) and the net
heat flux [W/m²] of a grey slab of thickness `L` in radiative equilibrium, with absorption
coefficient `κ`, scattering coefficient `σ_s` and Henyey–Greenstein asymmetry `g`. The walls
have emissivities `ε1`, `ε2`; a fraction `specularity` of the reflected power is mirrored,
the rest reflected diffusely. Returns `(T, q)`.
"""
function anisotropic_slab_equilibrium(κ::Real, σ_s::Real, g::Real, L::Real, T1::Real, T2::Real;
                                      ε1::Real = 1.0, ε2::Real = 1.0, specularity::Real = 0.0,
                                      Nx::Int = 32, n_half::Int = 32)
    β = κ + σ_s
    ω = σ_s / β
    μ, w = aslab_ordinates(n_half)
    n   = 2 * n_half
    K   = aslab_kernel(μ, w, Float64(ω), Float64(g), n_half - 1)
    Δτ  = β * L / Nx
    a   = exp.(-Δτ ./ abs.(μ))                       # transmission of one cell
    γ   = abs.(μ) ./ Δτ .* (1.0 .- a)                # cell average: Ī = (1 − γ) S + γ I_in
    pos = 1:n_half
    neg = n_half+1:n
    Ib1 = ASLAB_σ * T1^4 / π
    Ib2 = ASLAB_σ * T2^4 / π
    s   = Float64(specularity)

    # unknowns z = [cell-average intensities Ī (Nx × n); outgoing wall intensities I₁⁺, I₂⁻ (n_half each)]
    nI  = Nx * n
    nun = nI + n
    # one application of the affine map z ↦ Φ(z); `emit` switches the wall emission on (1) or off (0).
    # Returns the intensities arriving at the wall x = 0.
    function Φ!(out::Vector{Float64}, z::Vector{Float64}, emit::Float64)
        Ibar = reshape(view(z, 1:nI), Nx, n)
        new  = reshape(view(out, 1:nI), Nx, n)
        S    = Ibar * K'
        Iin  = z[nI .+ pos]
        for c in 1:Nx
            @views new[c, pos] .= (1.0 .- γ[pos]) .* S[c, pos] .+ γ[pos] .* Iin
            @views Iin .= a[pos] .* Iin .+ (1.0 .- a[pos]) .* S[c, pos]
        end
        arr2 = copy(Iin)                                       # arriving at x = L (μ > 0)
        Iin  = z[nI + n_half .+ pos]
        for c in Nx:-1:1
            @views new[c, neg] .= (1.0 .- γ[neg]) .* S[c, neg] .+ γ[neg] .* Iin
            @views Iin .= a[neg] .* Iin .+ (1.0 .- a[neg]) .* S[c, neg]
        end
        arr1 = copy(Iin)                                       # arriving at x = 0 (μ < 0), ordered as their mirrors
        q1 = 2.0 * sum(w[pos] .* μ[pos] .* arr1)               # incident flux / π at x = 0
        q2 = 2.0 * sum(w[pos] .* μ[pos] .* arr2)
        out[nI .+ pos]          .= emit * ε1 * Ib1 .+ (1.0 - ε1) .* ((1.0 - s) * q1 .+ s .* arr1)
        out[nI + n_half .+ pos] .= emit * ε2 * Ib2 .+ (1.0 - ε2) .* ((1.0 - s) * q2 .+ s .* arr2)
        return arr1
    end

    c0 = zeros(nun)
    Φ!(c0, zeros(nun), 1.0)
    A   = zeros(nun, nun)
    e   = zeros(nun)
    col = zeros(nun)
    for k in 1:nun
        e[k] = 1.0
        Φ!(col, e, 0.0)
        A[:, k] .= col
        e[k] = 0.0
    end
    z = (I - A) \ c0

    arr1 = Φ!(col, z, 1.0)
    Ibar = reshape(z[1:nI], Nx, n)
    T = (π .* 0.5 .* (Ibar * w) ./ ASLAB_σ) .^ 0.25            # σT⁴/π = ½ Σ w Ī
    q = 2π * (sum(w[pos] .* μ[pos] .* z[nI .+ pos]) - sum(w[pos] .* μ[pos] .* arr1))
    return T, q
end