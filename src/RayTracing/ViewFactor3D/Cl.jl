"""
    Cl(theta::Float64, almostZero::Float64)

Eq.(26) from the paper.
"""
function Cl(theta, almostZero)
    
    theta = mod(theta, 2*pi)
    chebArg = theta/pi - 1
    b = permutedims([1.865555351433979e-1, 6.269948963579612e-2, 3.139559104552675e-4,
            3.916780537368088e-6, 6.499672439854756e-8, 1.238143696612060e-9,
            5.586505893753557e-13])
    # Chebyshev polynomials of degrees 2*n+1 (n=1:6) found using the MATLAB sym command:
    # >> chebyshevT((2*(0:6)+1), sym(chebArg))
    T = permutedims([chebArg, 4*chebArg^3 - 3*chebArg,
            16*chebArg^5 - 20*chebArg^3 + 5*chebArg,
            64*chebArg^7 - 112*chebArg^5 + 56*chebArg^3 - 7*chebArg,
            256*chebArg^9 - 576*chebArg^7 + 432*chebArg^5 - 120*chebArg^3 + 9*chebArg,
            1024*chebArg^11 - 2816*chebArg^9 + 2816*chebArg^7 - 1232*chebArg^5 + 220*chebArg^3 - 11*chebArg,
            4096*chebArg^13 - 13312*chebArg^11 + 16640*chebArg^9 - 9984*chebArg^7 + 2912*chebArg^5 - 364*chebArg^3 + 13*chebArg])
    
    ClausenIntegral = ((theta - pi)*(2 + log((pi^2)/2)) + (2*pi - theta)*log((2*pi - theta)*(1 - almostZero) + almostZero) -
                        theta*log(theta*(1 - almostZero) + almostZero) + sum(b.*T))

    return ClausenIntegral

end