"""
    f(s::Float64, l::Float64, alpha::Float64,
        cosAlpha::Float64, sinAlpha::Float64,
        d::Float64, almostZero::Float64)

Eq.(22b) from the paper.
"""
function f(s, l, alpha, cosAlpha, sinAlpha, d, almostZero)

  s2 = s^2
  l2 = l^2
  d2 = d^2
  sinAlpha2 = sinAlpha^2
      
  wsqrt = sqrt(s2 + d2/sinAlpha2)
  psqrt = sqrt(l2 + d2/sinAlpha2)
  if abs(s + wsqrt) > 0
    wdim = s + wsqrt
  else 
    wdim = Float32(almostZero)
  end    
  if abs(l + psqrt) > 0
    pdim = l + psqrt
  else 
    pdim = Float32(almostZero)
  end

  F = ( (0.5*cosAlpha*(s2 + l2) - s*l)*log(s2 + l2 - 2*s*l*cosAlpha + d2) +
       s*sinAlpha*wsqrt*atan(sqrt(s2*sinAlpha2 + d2),(l - s*cosAlpha)) +
       l*sinAlpha*psqrt*atan(sqrt(l2*sinAlpha2 + d2),(s - l*cosAlpha)) + s*l +
       0.5*(d2/sinAlpha)*(imagLi_2((wdim/pdim), alpha, almostZero) +
       imagLi_2((pdim/wdim), alpha, almostZero) - 2*imagLi_2((wdim - 2*s)/pdim, (pi - alpha), almostZero))  )

    return Float32(F)
end