"""
    fParallel(s::Float64, l::Float64,
            d::Float64, almostZero::Float64)

Eq.(23) from the paper
"""
function fParallel(s::Float64, l::Float64,
                d::Float64, almostZero::Float64)
    
    if d == 0
        d = almostZero
    end
    
    sMinusl = s - l
    sMinusl2 = sMinusl^2
    s2 = s^2
    l2 = l^2
    d2 = d^2
    
    F = ( 0.5*(sMinusl2 - d2)*log(sMinusl2 + d2)
      - 2*sMinusl*d*acos(sMinusl/sqrt(s2 + l2 - 2*s*l + d2)) + s*l )

      return F
end