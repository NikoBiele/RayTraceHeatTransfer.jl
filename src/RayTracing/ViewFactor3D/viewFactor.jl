"""
    viewFactor(POLY_A::Matrix{Float64}, POLY_B::Matrix{Float64})

Original authors: Jacob A. Kerkhoff and Michael J. Wagner (MATLAB-implementation)
University of Wisconsin-Madison, Energy Systems Optimization Lab
https://github.com/uw-esolab/docs/tree/main/tools/viewfactor
Originally written for the paper:
"A Flexible Thermal Model for Solar Cavity Receivers Using Analytical View Factors"
https://doi.org/10.1115/ES2021-63810

Translated to the Julia Programming Language by: Nikolaj Maack Bielefeld
during studies at Aalborg University, Thermal Energy and Process Engineering.
The ray tracing part of the original function was not translated.
Published in RayTraceHeatTransfer.jl with permission from the original authors.

Analytical solution derivation:
Narayanaswamy, Arvind. "An analytic expression for radiation view 
factor between two arbitrarily oriented planar polygons." International
Journal of Heat and Mass Transfer 91 (2015): 841-847.

viewFactor(POLY_A, POLY_B) analytically computes the view factor from
POLY_A to POLY_B and from POLY_B to POLY_A. Input arguments are in
the form of 3x2 or Nx3 arrays, where each row corresponds to a vertex
of the polygon, the columns refer to the X,Y,Z coordinates of the vertices.
If Z coordinates are omitted, they are assumed to be zero.

The following must be true for both polygons:
1) polygons are planar (all vertices lie in the same plane)
2) polygons are simple (no self-intersecting polygons)
3) polygons are convex (in theory, concave polygons should work, but
    this remains untested)
"""
function viewFactor(POLY_A::Matrix{G}, POLY_B::Matrix{G}) where {G}

    # INITIALIZATION

    almostZero = 10*eps(G) # eps(Float32) # 1e-7
    halfTol = almostZero*10

    # CONFIRM CORRECT INPUTS AND FIND AREA: polygon A

    # also calculates unit normal vector to polygon
    dim_A = size(POLY_A)
    N = dim_A[1]  # number of vertices in polygon A
    # test for cartesian coordinates
    if dim_A[2] == 2
        POLY_A = [POLY_A, zeros(N,1)]
    elseif dim_A[2] != 3
        error("Input 1: Matrix dimensions incorrect.")
    end
    # test for coplanar vertices
    if N == 3    # test satisfied automatically
        n_A = cross((POLY_A[2,:] - POLY_A[1,:]), (POLY_A[3,:] - POLY_A[1,:]))
        nHat_A = n_A/norm(n_A)
        area_A = norm(n_A)/2   # area for triangles
    elseif N < 3 # not a polygon
        error("Input 1: Not enough coordinates to specify a polygon.")
    else
        n_A = cross((POLY_A[2,:] - POLY_A[1,:]), (POLY_A[3,:] - POLY_A[1,:]))
        nHat_A = n_A/norm(n_A)
        for i=4:N
            # the triple product of any combination of vertices must be zero
            #  for the polygon to be planar
            volume = abs(dot(n_A, (POLY_A[i,:] - POLY_A[1,:])))
            if volume > almostZero
                error("Input 1: Polygon vertices are not coplanar.")
            end
        end
        if N == 4   # area for arbitrary quadrilateral
            area_A = norm(cross((POLY_A[3,:] - POLY_A[1,:]), (POLY_A[4,:] - POLY_A[2,:])))/2
        else        # area for higher order polygons
            # http://geomalgorithms.com/a01-_area.html
            POLY_A_looped = [POLY_A ; permutedims(POLY_A[1,:])]
            toSum = [0; 0; 0]
            for i=1:N
                toSum = toSum + cross(POLY_A_looped[i,:], POLY_A_looped[i+1,:])
            end
            area_A = dot(nHat_A, toSum)/2
        end
    end

    # CONFIRM CORRECT INPUTS AND FIND AREA: polygon B

    # also calculates unit normal vector to polygon
    dim_B = size(POLY_B)
    M = dim_B[1]  # number of vertices in polygon B
    # test for cartesian coordinates
    if dim_B[2] == 2
        POLY_B = [POLY_B, zeros(M,1)]
    elseif dim_B[2] != 3
        error("Input 2: Matrix dimensions incorrect.")
    end
    # test for coplanar vertices
    if M == 3    # test satisfied automatically
        n_B = cross((POLY_B[2,:] - POLY_B[1,:]), (POLY_B[3,:] - POLY_B[1,:]))
        nHat_B = n_B/norm(n_B)
        area_B = norm(n_B)/2   # area for triangles
    elseif M < 3 # not a polygon
        error("Input 2: Not enough coordinates to specify a polygon.")
    else
        n_B = cross((POLY_B[2,:] - POLY_B[1,:]), (POLY_B[3,:] - POLY_B[1,:]))
        nHat_B = n_B/norm(n_B)
        for i=4:M
            # the triple product of any combination of vertices must be zero
            #  for the polygon to be planar
            volume = abs(dot(n_B, (POLY_B[i,:] - POLY_B[1,:])))
            if volume > almostZero
                error("Input 2: Polygon vertices are not coplanar.")
            end
        end
        if M == 4   # area for arbitrary quadrilateral
            area_B = norm(cross((POLY_B[3,:] - POLY_B[1,:]), (POLY_B[4,:] - POLY_B[2,:])))/2
        else        # area for higher order polygons
            # http://geomalgorithms.com/a01-_area.html
            POLY_B_looped = [POLY_B ; permutedims(POLY_B[1,:])]
            toSum = [0; 0; 0]
            for i=1:M
                toSum = toSum + cross(POLY_B_looped[i,:], POLY_B_looped[i+1,:])
            end
            area_B = dot(nHat_B, toSum)/2
        end
    end  

    # VIEW FACTOR ANALYICAL CALCULATION
        
    sumTerms = zeros(N,M)  # terms to sum to yield conductance
    skewPairs = zeros(N,M) # tracks which terms come from parallel edges (for debugging)
    for p = 1:M # loop through vertices of polygon B
        for i = 1:N  # loop through vertices of polygon A
            r_i = permutedims(POLY_A[i,:])
            r_p = permutedims(POLY_B[p,:])

            # loop pairings of vertices to cycle through edges
            if i < N
                r_j = permutedims(POLY_A[i+1,:])
            else # loop
                r_j = permutedims(POLY_A[1,:])
            end
            if p < M
                r_q = permutedims(POLY_B[p+1,:])
            else # loop
                r_q = permutedims(POLY_B[1,:])
            end

            # check for coincident vertices - nudge polygon B vertices if found
            if norm(r_i - r_p) < halfTol || norm(r_j - r_p) < halfTol
                r_p .+= almostZero
            elseif norm(r_i - r_q) < halfTol || norm(r_j - r_q) < halfTol
                r_q .+= almostZero
            end

            # determine parameterized coordinates for each edge, and minimum
            #  distance between edge rays (edges extended infinitely into space)
            dMin, sOrigin, sHat, lHat, lOrigin, skew = edgePairParameters(r_i, r_j, r_p, r_q, almostZero)

            if skew # if these edges are NOT parallel...
                # locate each vertex in the parameterized coordinate system
                s_i = dot((r_i - sOrigin), sHat)
                s_j = dot((r_j - sOrigin), sHat)
                l_p = dot((r_p - lOrigin), lHat)
                l_q = dot((r_q - lOrigin), lHat)

                skewPairs[i,p] = 1
                cosAlpha = dot(sHat, lHat) > 0.999 ? 0.999 : dot(sHat, lHat) < -0.999 ? -0.999 : dot(sHat, lHat)
                alpha = acos(cosAlpha)
                sinAlpha = sin(alpha)

                # Eq.(22a) from paper - calculate final terms that yield the 
                #  view factor when summed and divided by (4*pi*area)
                sumTerms[i,p] = (cosAlpha*(f(s_j, l_q, alpha, cosAlpha, sinAlpha, dMin, almostZero)
                    - f(s_i, l_q, alpha, cosAlpha, sinAlpha, dMin, almostZero)
                    - f(s_j, l_p, alpha, cosAlpha, sinAlpha, dMin, almostZero)
                    + f(s_i, l_p, alpha, cosAlpha, sinAlpha, dMin, almostZero)))
            else     # alternate expression for when alpha approaches zero  
                lHat = sHat # this is important for the parallel case
                # locate each vertex in the parameterized coordinate system
                s_i = dot((r_i - sOrigin), sHat)
                s_j = dot((r_j - sOrigin), sHat)
                l_p = dot((r_p - lOrigin), lHat)
                l_q = dot((r_q - lOrigin), lHat)

                skewPairs[i,p] = 0
                sumTerms[i,p] = dot(sHat, lHat)*(fParallel(s_j, l_q, dMin, almostZero)
                    - fParallel(s_i, l_q, dMin, almostZero) - fParallel(s_j, l_p, dMin, almostZero)
                    + fParallel(s_i, l_p, dMin, almostZero))
            end
        end
    end

    # "radiation conductance" : radUA = area_A*F_AB = area_B*F_BA
    radUA = abs(sum(sum(sumTerms)))/(4*pi)  

    F_AB = radUA/area_A
    F_BA = radUA/area_B

    return F_AB, F_BA, area_A, area_B
end