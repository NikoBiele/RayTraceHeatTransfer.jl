"""
    meshGeometry(subs::Vector{SubEnclosure},Ny::Int64,Nx::Int64) 

This function divides the geometry (any number of SubEnclosures)
into a 2D computational mesh over which the ray tracing is performed.
"""
function meshGeometry(subs::Vector{SubEnclosure},Ny::Int64,Nx::Int64) 
                    
    # The number of sub-enclosures.
    N_subs = length(subs)

    # First generate information describing the overall sub-enclosures.

    # Initialize arrays for each sub-enclosure.
    xs = Array{SVector{5,Float64}}(undef, N_subs) # Collection of the x-coordinates of each enclosure.
    ys = Array{SVector{5,Float64}}(undef, N_subs) # Collection of the y-coordinates of each enclosure.
    
    # Define boundary points of each sub-enclosure.
    for m = 1:N_subs
        # Define points in enclosure.
        Apoint = subs[m].point1
        Bpoint = subs[m].point2
        Cpoint = subs[m].point3
        Dpoint = subs[m].point4
        # x- and y-coordinates of this enclosure.
        xs_now = SVector(Apoint[1], Bpoint[1], Cpoint[1], Dpoint[1], Apoint[1])
        ys_now = SVector(Apoint[2], Bpoint[2], Cpoint[2], Dpoint[2], Apoint[2])
        xs[m] = xs_now # Save the vector of points.
        ys[m] = ys_now
    end

    # Split each sub-enclosure into a specified number cells.

    # Initialize arrays to hold information about each cell.
    xPoints = zeros(Nx+1,Ny+1)
    yPoints = zeros(Nx+1,Ny+1)
    point1 = Matrix{SVector{2,Float64}}(undef, Nx, N_subs*Ny)
    point2 = Matrix{SVector{2,Float64}}(undef, Nx, N_subs*Ny)
    point3 = Matrix{SVector{2,Float64}}(undef, Nx, N_subs*Ny)
    point4 = Matrix{SVector{2,Float64}}(undef, Nx, N_subs*Ny)

    # Loop over all sub-enclosures.
    for k = 1:N_subs

        xs_now = xs[k] # x- and y-points of the boundary of this SubEnclosure
        ys_now = ys[k]

        # Improved algorithm!

        # Get change of coordinates at outer dimensions.
        # x first.
        deltaXbot = xs_now[2]-xs_now[1]
        deltaXtop = xs_now[4]-xs_now[3]
        deltaXleft = xs_now[5]-xs_now[4]
        # then y.
        deltaYBot = ys_now[1]-ys_now[2]
        deltaYRight = ys_now[2] - ys_now[3]
        deltaYLeft = ys_now[4] - ys_now[1]

        # Loop to create sub-divisions (cells) of each sub-enclosure.
        for m = 1:Ny+1 # loop over y

            # We move the reference point in each iteration.
            refmoveXleft = (m-1)*deltaXleft/Ny
            refmoveXright = deltaXbot - (m-1)*(deltaXbot+deltaXtop)/Ny

            # Start at the bottom left point in any enclosure.
            for n = 1:Nx+1

                # Change in x and y.
                # Add or subtract a bit of y.
                refmoveYdown = (n-1)*deltaYBot/Nx
                refmoveYup = deltaYLeft - (n-1)*(deltaYLeft+deltaYRight)/Nx

                # First term = bottom left point.
                # Second term = how much we move the reference.
                # Third term = how much we move right in each iteration.
                xPoints[n,m] = xs_now[1] - refmoveXleft + (n-1)*(refmoveXright)/Nx
                yPoints[n,m] = ys_now[1] - refmoveYdown + (m-1)*(refmoveYup)/Ny
            end
        end
        
        # Put the points into the matrices to work with ray tracing.
        for m = 1:Ny
            for n = 1:Nx
                point1[n,m+(k-1)*Ny] = SVector(xPoints[n,m], yPoints[n,m])
                point2[n,m+(k-1)*Ny] = SVector(xPoints[n+1,m], yPoints[n+1,m])
                point3[n,m+(k-1)*Ny] = SVector(xPoints[n+1,m+1], yPoints[n+1,m+1])
                point4[n,m+(k-1)*Ny] = SVector(xPoints[n,m+1], yPoints[n,m+1])
            end
        end
    end

    # miscellaneous useful numbers related to geometry
    N_surfs = 2*Nx+2*Ny*N_subs
    N_vols = Nx*Ny*N_subs

    return point1, point2, point3, point4, N_surfs, N_vols

end