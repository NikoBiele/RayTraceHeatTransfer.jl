"""
    Cl3D(theta)

Clausen function Cl₂(θ) = −∫₀^θ log|2 sin(t/2)| dt, Eq. (26) of the paper.

The paper gives a 7-term Chebyshev fit for Eq. (26), which is accurate to only
about 3e-11 even in exact arithmetic. This replaces it by the series

    Cl₂(θ) = θ − θ log|θ| + Σₖ |B₂ₖ| θ^(2k+1) / (2k (2k+1) (2k)!),   |θ| ≤ π,

after reducing θ to (−π, π] (Cl₂ is odd and 2π-periodic). Successive terms
fall by about (θ/2π)² ≤ 1/4, and the sum stops once they drop below eps(T),
so the result is accurate to full precision in any floating-point type T
(about 27 terms in Float64, 130 at 256-bit BigFloat). For BigFloat, set the
precision globally with `setprecision` before calling.
"""
function Cl3D(theta)
    T = float(typeof(theta))
    θ = mod(T(theta), 2 * T(pi))            # reduce to [0, 2π)
    θ > T(pi) && (θ -= 2 * T(pi))           # then to (−π, π]
    iszero(θ) && return zero(T)             # Cl₂(0) = 0, since θ log|θ| → 0

    a   = _clausen_coeffs(T)                # series coefficients, rounded to T
    x   = θ * θ                             # θ²
    p   = θ * x                             # θ^(2k+1), starting at k = 1
    tol = eps(T) / 8                        # stop once a term no longer counts
    acc = θ - θ * log(abs(θ))               # the closed-form part
    for k in eachindex(a)
        term = a[k] * p                     # k-th series term
        acc += term
        abs(term) < tol && break            # terms decrease monotonically
        p *= x                              # next odd power of θ
    end
    return acc
end

# Tangent numbers T₁, …, T_K as exact integers (Brent & Harvey 2011).
# They give the Bernoulli numbers via |B₂ₖ| = 2k Tₖ / (4ᵏ (4ᵏ − 1)).
function _tangent_numbers(K::Int)
    t = zeros(BigInt, K)
    t[1] = 1
    for k in 2:K
        t[k] = (k - 1) * t[k-1]
    end
    for k in 2:K, j in k:K
        t[j] = (j - k) * t[j-1] + (j - k + 2) * t[j]
    end
    return t
end

# Series coefficients aₖ = |B₂ₖ| / (2k (2k+1) (2k)!) = Tₖ / (4ᵏ (4ᵏ − 1) (2k+1) (2k)!),
# formed from exact integers and rounded once to T. At |θ| = π the terms fall
# by about 1/4 each, so K = bits/2 terms reach eps(T); 8 more are a margin.
function _clausen_coeffs_exact(::Type{T}, bits::Int) where {T}
    K = cld(bits, 2) + 8
    t = _tangent_numbers(K)
    a = Vector{T}(undef, K)
    for k in 1:K
        den  = big(4)^k * (big(4)^k - 1) * (2k + 1) * factorial(big(2k))  # exact denominator
        a[k] = T(BigFloat(t[k]) / BigFloat(den))                          # one rounding to T
    end
    return a
end

# Float64 table built once at load, so threaded callers never take the lock.
const _CLAUSEN_F64 = _clausen_coeffs_exact(Float64, 53)

# Other types (BigFloat, Double64, …): built on first use, cached per type and precision.
const _CLAUSEN_CACHE = Dict{Tuple{DataType,Int},Any}()
const _CLAUSEN_LOCK  = ReentrantLock()

_clausen_coeffs(::Type{Float64}) = _CLAUSEN_F64

function _clausen_coeffs(::Type{T}) where {T}
    bits = ceil(Int, -log2(eps(one(T))))    # working precision of T in bits
    coeffs = lock(_CLAUSEN_LOCK) do
        get!(() -> _clausen_coeffs_exact(T, bits), _CLAUSEN_CACHE, (T, bits))
    end
    return coeffs::Vector{T}
end