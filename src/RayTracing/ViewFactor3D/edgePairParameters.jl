"""
    edgePairParameters(Po::Matrix{Float64}, Pf::Matrix{Float64},
                        Qo::Matrix{Float64}, Qf::Matrix{Float64},
                        almostZero::Float64)

http://geomalgorithms.com/a07-_distance.html
find shortest distance D between line Po+s*u and Qo+t*v for initial 
points Po and Qo, parameters s and t, and vectors u and v.
"""
function edgePairParameters(Po, Pf, Qo, Qf,almostZero)
    
    u = Pf - Po
    v = Qf - Qo
    w = Po - Qo
    
    Plength = norm(u)
    Qlength = norm(v)
    u = u/Plength  # make these unit vectors
    v = v/Qlength
    
    a = 1 # dot(u, u)
    b = dot(u, v)
    c = 1 # dot(v, v)
    d = dot(u, w)
    e = dot(v, w)
    
    den = a*c - b^2
    
    # calculate shortest distance between edge rays
    if den > almostZero
        skew = true
        s = (b*e - c*d)/den
        l = (a*e - b*d)/den
        D = norm(w + s*u - l*v)
    else # origin is arbitrary if lines are parallel
        skew = false
    #     s = 1.5*Plength
    #     l = 1.5*Qlength
        s = 0
        l = e/c
        D = norm(w - (e/c)*v)
    end
    
    # see Fig 5 in this paper:
    #   Narayanaswamy, Arvind. "An analytic expression for radiation view 
    #   factor between two arbitrarily oriented planar polygons." International
    #   Journal of Heat and Mass Transfer 91 (2015): 841-847.
    # for description of why these values are calculated in this way.
    
    # parameter origin is location on edge ray where distance between edges has
    #  its smallest value 
    sOrigin = Po + u*s
    lOrigin = Qo + v*l
    
    s_toEnd = norm(Pf - sOrigin)
    l_toEnd = norm(Qf - lOrigin)
    
    # unit vectors point from parameter origin to furthest of the two vertices
    if abs(s) < s_toEnd
        sHat = (Pf - sOrigin)/norm(Pf - sOrigin)
    else
        sHat = (Po - sOrigin)/norm(Po - sOrigin)
    end
    if abs(l) < l_toEnd
        lHat = (Qf - lOrigin)/norm(Qf - lOrigin)
    else
        lHat = (Qo - lOrigin)/norm(Qo - lOrigin)
    end
    
    return D, sOrigin, sHat, lHat, lOrigin, skew
end