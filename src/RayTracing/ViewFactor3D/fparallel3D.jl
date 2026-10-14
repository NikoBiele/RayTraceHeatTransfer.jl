"""
    fParallel3D(s, l, d, almostZero)

Eq. (23) from the paper: the edge-pair term for parallel edges a distance d apart.

Written exactly: acos((s − l)/√((s − l)² + d²)) is evaluated as atan(d, s − l),
which holds for d ≥ 0 and needs neither clamping nor a regularised square root,
and the logarithmic term takes its limit 0 where (s − l)² + d² = 0. The result
is accurate to the precision of the input type. `almostZero` is no longer used
and is kept only so the callers are unchanged.
"""
function fParallel3D(s, l, d, almostZero)
    a  = s - l                                     # offset along the common edge direction
    a2 = a^2                                       # (s − l)², formed without cancellation
    d2 = d^2
    r2 = a2 + d2                                   # squared distance between the edge points
    logTerm = iszero(r2) ? zero(r2) : (a2 - d2) * log(r2) / 2   # limit 0 as r² → 0
    return logTerm - 2 * a * d * atan(d, a) + s * l             # acos(a/√r²) = atan(d, a), d ≥ 0
end