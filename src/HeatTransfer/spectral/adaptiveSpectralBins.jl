# Adaptive spectral binning by κ-level sets. Given absorption-coefficient
# samples κ(λ), bins are intervals in κ, refined by greedy bisection until
# the energy-weighted first-pass transmission error
#
#     E(L; T) = Σ_bins f_bin(T) · | ⟨exp(−κL)⟩_bin,T − exp(−κ_ref L) |
#
# is below `tol` for every L in the path-length range and every T in the
# temperature range. Bin membership is decided per sample (the interpolant
# only locates the crossing wavelengths between samples), so a poorly chosen
# kernel can misplace an edge but never change which samples form a bin.

struct _KappaBin
    lo::Float64
    hi::Float64
    samples::Vector{Int}
    κ_ref::Float64
    err::Matrix{Float64}          # length(Ls) × length(temps)
end

_planck_spectral(λ::Real, T::Real) = C1 / (λ^5 * expm1(C2 / (λ * T)))

function _trapezoid_weights(λ::AbstractVector)
    n = length(λ)
    w = zeros(n)
    w[1] = (λ[2] - λ[1]) / 2
    w[n] = (λ[n] - λ[n-1]) / 2
    @inbounds for i in 2:n-1
        w[i] = (λ[i+1] - λ[i-1]) / 2
    end
    return w
end

function _error_matrix(samples, κ_ref, κ, Bw, Btot, Ls)
    err = zeros(length(Ls), size(Bw, 2))
    isempty(samples) && return err
    κs = κ[samples]
    for t in 1:size(Bw, 2)
        bw = Bw[samples, t]
        Wb = sum(bw)
        Wb == 0 && continue
        for (m, L) in enumerate(Ls)
            T_exact = sum(bw .* exp.(-κs .* L)) / Wb
            err[m, t] = (Wb / Btot[t]) * abs(T_exact - exp(-κ_ref * L))
        end
    end
    return err
end

function _make_kbin(lo, hi, samples, κ, Bw, Btot, Ls)
    κ_ref = isempty(samples) ? NaN : exp(sum(log.(κ[samples])) / length(samples))
    return _KappaBin(lo, hi, samples, κ_ref, _error_matrix(samples, κ_ref, κ, Bw, Btot, Ls))
end

# Wavelength where the interpolant of log10 κ over log10 λ crosses log10 ℓ
# between samples i and i+1: bisection on the interpolant when it brackets
# the level, log-log linear interpolation between the samples otherwise.
function _crossing(u, v, itp, i, ℓ)
    target = log10(ℓ)
    a, b = u[i], u[i+1]
    fa, fb = itp(a) - target, itp(b) - target
    if fa * fb < 0
        for _ in 1:60
            c = (a + b) / 2
            fc = itp(c) - target
            if fa * fc <= 0
                b, fb = c, fc
            else
                a, fa = c, fc
            end
        end
        return 10.0^((a + b) / 2)
    end
    t = (target - v[i]) / (v[i+1] - v[i])
    return 10.0^(u[i] + clamp(t, 0.0, 1.0) * (u[i+1] - u[i]))
end

"""
    adaptiveSpectralBins(λ, κ; tol, L_range, T_range,
                         τ_window = (0.01, 100.0), scale_range = (1.0, 1.0),
                         r0 = 2.0, max_bins = 5000, kernel = :b7) -> PiecewiseBands

Build κ-level-set spectral bins from absorption-coefficient samples `κ`
[1/m] at wavelengths `λ` [m] (increasing), refined until the energy-weighted
first-pass transmission error is below `tol` for all path lengths in
`L_range` (cell size .. enclosure size, m) and temperatures in `T_range` (K).

Per-element spectra of the form κ_e(λ) = s_e·κ(λ) are supported by giving
the range of scale factors in `scale_range`; bins are then valid for all
elements. `τ_window` sets the initial coarse levels (ratio `r0`) in optical
depth; refinement is driven by the tolerance alone. `kernel` is the
ConvolutionInterpolations kernel used to locate crossing wavelengths between
samples.

The returned model carries `κ_ref` (use `face.kappa_g = model.κ_ref .* s_e`)
and `achieved_error`, the bound actually reached.
"""
function adaptiveSpectralBins(λ::AbstractVector, κ::AbstractVector;
                              tol::Real, L_range::Tuple{<:Real,<:Real}, T_range::Tuple{<:Real,<:Real},
                              τ_window::Tuple{<:Real,<:Real} = (0.01, 100.0),
                              scale_range::Tuple{<:Real,<:Real} = (1.0, 1.0),
                              r0::Real = 2.0, max_bins::Int = 5000, kernel::Symbol = :b7)
    λ = Float64.(λ); κ = Float64.(κ)
    n = length(λ)
    n == length(κ) || throw(ArgumentError("λ and κ must have the same length"))
    n >= 4 || throw(ArgumentError("need at least 4 samples"))
    all(diff(λ) .> 0) || throw(ArgumentError("λ must be strictly increasing"))
    all(λ .> 0) || throw(ArgumentError("λ must be positive"))
    all(κ .> 0) || throw(ArgumentError("κ must be positive (use a small floor for transparent regions)"))
    tol > 0 || throw(ArgumentError("tol must be positive"))

    L_min = Float64(L_range[1]) * Float64(scale_range[1])
    L_max = Float64(L_range[2]) * Float64(scale_range[2])
    0 < L_min <= L_max || throw(ArgumentError("invalid L_range / scale_range"))
    Ls = 10 .^ range(log10(L_min), log10(L_max), length = 31)
    temps = [Float64(T_range[1]), Float64(T_range[2])]
    all(temps .> 0) || throw(ArgumentError("T_range must be positive"))

    # Planck × quadrature weights at both temperatures
    w = _trapezoid_weights(λ)
    Bw = hcat([_planck_spectral.(λ, T) .* w for T in temps]...)
    Btot = vec(sum(Bw, dims = 1))

    # ---- initial κ-intervals: coarse window levels clipped to the data range
    κmin, κmax = extrema(κ)
    κ_top = κmax * (1 + 1e-12)
    κ_lo = Float64(τ_window[1]) / L_max
    κ_hi = Float64(τ_window[2]) / L_min
    J = max(0, ceil(Int, log(κ_hi / κ_lo) / log(r0)))
    levels = [κ_lo * r0^j for j in 0:J]
    edges = sort(unique(vcat(κmin, filter(ℓ -> κmin < ℓ < κmax, levels), κ_top)))

    b_of = [searchsortedlast(edges, k) for k in κ]
    bins = _KappaBin[]
    for b in 1:length(edges)-1
        push!(bins, _make_kbin(edges[b], edges[b+1], findall(==(b), b_of), κ, Bw, Btot, Ls))
    end
    total = zeros(length(Ls), length(temps))
    for b in bins
        total .+= b.err
    end

    # ---- greedy bisection of the bin with the largest peak contribution
    while maximum(total) > tol
        length(bins) >= max_bins &&
            error("adaptiveSpectralBins: $(max_bins) bins reached with error $(maximum(total)) > tol; raise max_bins, tol, or supply better-resolved κ samples")
        ib = argmax([maximum(b.err) for b in bins])
        b = bins[ib]
        length(b.samples) < 2 &&
            error("adaptiveSpectralBins: cannot reach tol = $(tol) (best $(maximum(total))): a bin is down to one sample; the κ samples do not resolve the spectrum finely enough")
        mid = sqrt(b.lo * b.hi)
        s_lo = filter(i -> κ[i] < mid, b.samples)
        if isempty(s_lo) || length(s_lo) == length(b.samples)
            logs = sort(log.(κ[b.samples]))
            mid = exp(logs[cld(length(logs), 2)])
            s_lo = filter(i -> κ[i] < mid, b.samples)
            (isempty(s_lo) || length(s_lo) == length(b.samples)) &&
                error("adaptiveSpectralBins: cannot split a bin whose samples share one κ value; tol = $(tol) unreachable (best $(maximum(total)))")
        end
        s_hi = filter(i -> κ[i] >= mid, b.samples)
        new_lo = _make_kbin(b.lo, mid, s_lo, κ, Bw, Btot, Ls)
        new_hi = _make_kbin(mid, b.hi, s_hi, κ, Bw, Btot, Ls)
        total .-= b.err
        total .+= new_lo.err
        total .+= new_hi.err
        bins[ib] = new_lo
        push!(bins, new_hi)
    end
    achieved = maximum(total)

    # ---- bins in κ order; sample → bin index
    sort!(bins, by = b -> b.lo)
    bins = filter(b -> !isempty(b.samples), bins)
    K = length(bins)
    bin_of = zeros(Int, n)
    for (k, b) in enumerate(bins)
        bin_of[b.samples] .= k
    end

    # ---- pieces: crossing wavelengths where consecutive samples change bin
    u = log10.(λ)
    v = log10.(κ)
    itp = convolution_interpolation(u, v; kernel = kernel)
    piece_edges = [λ[1]]
    piece_bin = Int[]
    cur = bin_of[1]
    for i in 1:n-1
        nxt = bin_of[i+1]
        nxt == cur && continue
        # bins strictly between cur and nxt are crossed in order; each crossing
        # is at the shared κ-edge of two consecutive bins
        step = nxt > cur ? 1 : -1
        k = cur
        while k != nxt
            ℓ = step == 1 ? bins[k].hi : bins[k].lo
            λc = clamp(_crossing(u, v, itp, i, ℓ), λ[i], λ[i+1])
            λc = max(λc, piece_edges[end] * (1 + 1e-12))     # keep edges strictly increasing
            push!(piece_edges, λc)
            push!(piece_bin, k)
            k += step
        end
        cur = nxt
    end
    push!(piece_edges, λ[n])
    push!(piece_bin, cur)

    return PiecewiseBands(piece_edges, piece_bin,
                          [b.κ_ref for b in bins], [b.lo for b in bins], [b.hi for b in bins],
                          achieved)
end