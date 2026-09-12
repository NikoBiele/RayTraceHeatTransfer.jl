# Spectral emission-fraction models.
#
# A spectral domain carries one `spectral_model`, which decides how an
# element's blackbody emission at temperature T is divided over the
# n_spectral_bins bins (and how that division changes with T, for the
# Newton Jacobian). Each concrete model implements:
#
#     n_bins(m)                              number of bins
#     validate(m, n_bins_domain)             consistency checks, error on failure
#     describe(m)                            one-line text for show()
#     fill_bin_fractions!(row, m, T)         row[k] = f_k(T), sums to 1
#     bin_fraction_derivatives!(row, m, T)   row[k] = ∂f_k/∂T
#
# Everything else in the spectral solvers is model-agnostic.

"""
    PlanckBands(limits)

Contiguous wavelength bands: bin k spans [limits[k], limits[k+1]] (in m), the
first band extended down to λ = 0 and the last up to λ = ∞, so the Planck
fractions sum to 1 at every temperature. `limits` must be strictly increasing
and positive, with at least 4 entries (3 bins).
"""
struct PlanckBands <: AbstractSpectralModel
    limits::Vector{Float64}
end

"""
    ConstantWeights(a)

Temperature-independent bin weights. `a` must be non-negative and sum to 1.
"""
struct ConstantWeights <: AbstractSpectralModel
    a::Vector{Float64}
end

# ---------------------------------------------------------------- n_bins ---
n_bins(m::PlanckBands)     = length(m.limits) - 1
n_bins(m::ConstantWeights) = length(m.a)

# -------------------------------------------------------------- validate ---
function validate(m::PlanckBands, n_bins_domain::Int)
    lims = m.limits
    length(lims) >= 4 || error("PlanckBands: limits must have at least 4 values (defining 3 bins)")
    all(lims .> 0) || error("PlanckBands: limits must all be positive (wavelengths > 0)")
    all(diff(lims) .> 0) || error("PlanckBands: limits must be strictly increasing (no duplicates)")
    n_bins(m) == n_bins_domain ||
        error("PlanckBands defines $(n_bins(m)) bins but the domain has $(n_bins_domain) spectral bins")
    return nothing
end

function validate(m::ConstantWeights, n_bins_domain::Int)
    a = m.a
    length(a) == n_bins_domain ||
        error("ConstantWeights has $(length(a)) entries but the domain has $(n_bins_domain) spectral bins")
    all(a .>= 0) || error("ConstantWeights: weights must be non-negative")
    isapprox(sum(a), 1.0; atol = 1e-8) ||
        error("ConstantWeights: weights must sum to 1 (sum = $(sum(a)))")
    return nothing
end

# -------------------------------------------------------------- describe ---
describe(m::PlanckBands) =
    string("λ ∈ [", round(first(m.limits); sigdigits = 3), ", ",
           round(last(m.limits); sigdigits = 3), "] (Planck bands)")
describe(m::ConstantWeights) =
    string("constant weights ", round.(m.a; sigdigits = 3))

# ---------------------------------------------------- emission fractions ---
"""
    fill_bin_fractions!(row, m::PlanckBands, T)

Bin k gets F(limits[k+1]·T) − F(limits[k]·T) with F the cumulative blackbody
fraction; first bin from F = 0, last bin up to F = 1.
"""
function fill_bin_fractions!(row::AbstractVector, m::PlanckBands, T::Real)
    K = n_bins(m)
    F_prev = 0.0
    for k in 1:K
        if k == K
            row[k] = 1.0 - F_prev
        else
            F_curr = emitFracBlackBodySpectrum(m.limits, T, k + 1)
            row[k] = F_curr - F_prev
            F_prev = F_curr
        end
    end
    return nothing
end

function fill_bin_fractions!(row::AbstractVector, m::ConstantWeights, T::Real)
    row .= m.a
    return nothing
end

# ------------------------------------------------- fraction derivatives ---
"""
    bin_fraction_derivatives!(row, m::PlanckBands, T)

∂f_k/∂T = λ_{k+1}·F′(λ_{k+1}T) − λ_k·F′(λ_k T), with F′(x) = dF/dx evaluated
at x = λT, and the first/last edges extended to 0/∞ (zero derivative).
"""
function bin_fraction_derivatives!(row::AbstractVector, m::PlanckBands, T::Real)
    K = n_bins(m)
    lam = m.limits
    if T <= 0
        row .= 0.0
        return nothing
    end
    for k in 1:K
        dF_lo = (k == 1) ? 0.0 : lam[k]   * dF_blackbody_dlambdaT(lam[k]   * T)
        dF_hi = (k == K) ? 0.0 : lam[k+1] * dF_blackbody_dlambdaT(lam[k+1] * T)
        row[k] = dF_hi - dF_lo
    end
    return nothing
end

function bin_fraction_derivatives!(row::AbstractVector, m::ConstantWeights, T::Real)
    row .= 0.0
    return nothing
end

# ---------------------------------------------------------------------------
# PiecewiseBands — bins that are unions of wavelength intervals
# ---------------------------------------------------------------------------

"""
    PiecewiseBands(edges, piece_bin, κ_ref, κ_lo, κ_hi, achieved_error)

Bins defined by κ-level sets: the wavelength axis is partitioned into pieces
[edges[p], edges[p+1]], p = 1..P, and piece p belongs to bin `piece_bin[p]`.
A bin's emission fraction at temperature T is the sum of the Planck fractions
of its pieces, with the first piece extended down to λ = 0 and the last up to
λ = ∞, so the fractions sum to 1 at every temperature. `κ_ref[k]` is the
reference absorption coefficient of bin k (used as `kappa_g[k]` up to a
per-element scale), `κ_lo[k] ≤ κ < κ_hi[k]` its κ-interval, and
`achieved_error` the transmission-error bound reported by the constructor
(`adaptiveSpectralBins`). Usually built by that constructor, not by hand.
"""
struct PiecewiseBands <: AbstractSpectralModel
    edges::Vector{Float64}       # P + 1 increasing wavelengths [m]
    piece_bin::Vector{Int}       # P bin indices
    κ_ref::Vector{Float64}       # K
    κ_lo::Vector{Float64}        # K
    κ_hi::Vector{Float64}        # K
    achieved_error::Float64
end

n_bins(m::PiecewiseBands) = length(m.κ_ref)
n_pieces(m::PiecewiseBands) = length(m.piece_bin)

function validate(m::PiecewiseBands, n_bins_domain::Int)
    K, P = n_bins(m), n_pieces(m)
    K == n_bins_domain ||
        error("PiecewiseBands defines $(K) bins but the domain has $(n_bins_domain) spectral bins")
    length(m.edges) == P + 1 || error("PiecewiseBands: edges must have one more entry than piece_bin")
    all(m.edges .> 0) || error("PiecewiseBands: edges must be positive wavelengths")
    all(diff(m.edges) .> 0) || error("PiecewiseBands: edges must be strictly increasing")
    all(1 .<= m.piece_bin .<= K) || error("PiecewiseBands: piece_bin entries must be in 1:$(K)")
    length(unique(m.piece_bin)) == K || error("PiecewiseBands: every bin must own at least one piece")
    (length(m.κ_lo) == K && length(m.κ_hi) == K) ||
        error("PiecewiseBands: κ_lo and κ_hi must have one entry per bin")
    all(m.κ_ref .>= 0) || error("PiecewiseBands: κ_ref must be non-negative")
    return nothing
end

describe(m::PiecewiseBands) =
    string(n_bins(m), " κ-level bins over ", n_pieces(m), " pieces, κ_ref ∈ [",
           round(minimum(m.κ_ref); sigdigits = 3), ", ",
           round(maximum(m.κ_ref); sigdigits = 3), "], bound ",
           round(m.achieved_error; sigdigits = 2))

# Function barriers: `itp` is fetched once from the untyped cache and passed
# in, so the loops specialise on its concrete type.
function _piecewise_fractions!(row::AbstractVector, m::PiecewiseBands, T::Real, itp)
    row .= 0.0
    P = n_pieces(m)
    F_prev = 0.0                                  # first piece extends to λ = 0
    @inbounds for p in 1:P
        F_next = (p == P) ? 1.0 : planck_F(m.edges[p+1] * T, itp)   # last piece extends to λ = ∞
        row[m.piece_bin[p]] += F_next - F_prev
        F_prev = F_next
    end
    return nothing
end

function _piecewise_derivatives!(row::AbstractVector, m::PiecewiseBands, T::Real)
    row .= 0.0
    P = n_pieces(m)
    d_prev = 0.0                                  # d/dT of F at λ = 0 is 0
    @inbounds for p in 1:P
        d_next = (p == P) ? 0.0 : m.edges[p+1] * dF_blackbody_dlambdaT(m.edges[p+1] * T)
        row[m.piece_bin[p]] += d_next - d_prev
        d_prev = d_next
    end
    return nothing
end

function fill_bin_fractions!(row::AbstractVector, m::PiecewiseBands, T::Real)
    if T <= 0
        row .= 0.0
        return nothing
    end
    _piecewise_fractions!(row, m, T, planck_table())
    return nothing
end

function bin_fraction_derivatives!(row::AbstractVector, m::PiecewiseBands, T::Real)
    if T <= 0
        row .= 0.0
        return nothing
    end
    _piecewise_derivatives!(row, m, T)
    return nothing
end