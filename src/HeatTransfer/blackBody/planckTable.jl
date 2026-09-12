# Tabulated cumulative blackbody fraction F(x), x = λT [m·K].
#
# The piecewise spectral model evaluates F at two edges per piece for every
# element at every temperature update — tens of thousands of evaluations per
# element for a line spectrum — so the series in emitFracBlackBodySpectrum is
# replaced by a table in u = log10(x), interpolated with a :b7 convolution
# kernel (7th-order accurate; 4001 nodes over five decades give errors far
# below 1e-10). Outside the table F is 0 below and 1 above, which is exact to
# better than 1e-8 at the chosen limits. The derivative dF/dx keeps its closed
# form (dF_blackbody_dlambdaT).
#
# The table is built on first use and cached; hot loops fetch it once and
# pass it through a function barrier so the interpolant type is concrete.

const _PLANCK_U_LO = -4.5          # x = 3.2e-5 m·K  → F ≈ 1e-190
const _PLANCK_U_HI =  0.5          # x = 3.2   m·K  → 1 − F ≈ 5e-9
const _PLANCK_X_LO = 10.0^_PLANCK_U_LO
const _PLANCK_X_HI = 10.0^_PLANCK_U_HI
const _PLANCK_TABLE = Ref{Any}(nothing)

"""
    planck_table()

The cached interpolant of F(x) in u = log10(x). Built on first call.
"""
function planck_table()
    t = _PLANCK_TABLE[]
    t === nothing || return t
    u = collect(range(_PLANCK_U_LO, _PLANCK_U_HI, length = 4001))
    F = [emitFracBlackBodySpectrum((10.0^ui,), 1.0, 1) for ui in u]
    itp = convolution_interpolation(u, F; kernel = :b7)
    _PLANCK_TABLE[] = itp
    return itp
end

"""
    planck_F(x, itp)

Cumulative blackbody fraction F(0→x) at x = λT [m·K], from the table `itp`
returned by `planck_table()`. Clamped to 0 / 1 outside the tabulated range.
"""
@inline function planck_F(x::Real, itp)
    x <= _PLANCK_X_LO && return 0.0
    x >= _PLANCK_X_HI && return 1.0
    return itp(log10(x))
end