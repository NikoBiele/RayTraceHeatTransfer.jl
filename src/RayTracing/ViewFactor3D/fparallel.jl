"""
    fParallel(s::Float64, l::Float64,
            d::Float64, almostZero::Float64)

Eq.(23) from the paper
"""
function fParallel(s, l, d, almostZero)
    
    if d == 0
        d = Float32(almostZero)
    end
    
    sMinusl = s - l
    sMinusl2 = sMinusl^2
    s2 = s^2
    l2 = l^2
    d2 = d^2
    term = (sMinusl/sqrt(s2 + l2 - 2*s*l + d2 + almostZero)) 
    term = term >= 0.999999 ? 0.999999 : term <= -0.999999 ? -0.999999 : term
    
    F = ( 0.5*(sMinusl2 - d2)*log(sMinusl2 + d2)
            - 2*sMinusl*d*acos(term) + s*l )

      return F
end