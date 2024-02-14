function geometry(yLayersHeight::Vector{Float64},xLayersWidth::Matrix{Float64},
                    Ny::Int64,Nx::Int64,displayGeometry::Bool)
    # This function divides the geometry into a 2D computational mesh
    # over which the ray tracing is performed.    
                    
    # function to plot geometry
    function trapezoid(point1, point2, point3, point4)
        Shape([point1[1], point2[1], point3[1], point4[1]],[point1[2], point2[2], point3[2], point4[2]])
    end

    # The number of sub-enclosures (height layers minus one).
    N_subs = length(yLayersHeight)-1

    ##### First generate information describing the overall sub-enclosures.

    # Initialize arrays for each sub-enclosure.
    xs = Array{SVector{5,Float64}}(undef, N_subs) # Collection of the x-coordinates of each enclosure.
    ys = Array{SVector{5,Float64}}(undef, N_subs) # Collection of the y-coordinates of each enclosure.
    A = Array{SVector{2,Float64}}(undef, N_subs) # Bottom left points.
    B = Array{SVector{2,Float64}}(undef, N_subs) # Bottom right points.
    C = Array{SVector{2,Float64}}(undef, N_subs) # Top right points.
    D = Array{SVector{2,Float64}}(undef, N_subs) # Top left points.
    
    # Define boundary points of each sub-enclosure.
    for m = 1:N_subs
        # Define points in enclosure.
        Apoint = SVector(xLayersWidth[1,m], yLayersHeight[m])
        Bpoint = SVector(xLayersWidth[2,m], yLayersHeight[m])
        Cpoint = SVector(xLayersWidth[2,m+1], yLayersHeight[m+1])
        Dpoint = SVector(xLayersWidth[1,m+1], yLayersHeight[m+1])
        # Save points.
        A[m] = Apoint
        B[m] = Bpoint
        C[m] = Cpoint
        B[m] = Dpoint
        # x- and y-coordinates of this enclosure.
        xs_now = SVector(Apoint[1], Bpoint[1], Cpoint[1], Dpoint[1], Apoint[1])
        ys_now = SVector(Apoint[2], Bpoint[2], Cpoint[2], Dpoint[2], Apoint[2])
        xs[m] = xs_now # Save the vector of points.
        ys[m] = ys_now
    end

    ##### Split each sub-enclosure into a specified number cells.

    # Initialize arrays to hold information about each cell.
    xPoints = zeros(Nx+1,Ny+1)
    yPoints = zeros(Nx+1,Ny+1)
    point1 = Array{SVector{2,Float64}}(undef, Nx, Ny*N_subs)
    point2 = Array{SVector{2,Float64}}(undef, Nx, Ny*N_subs)
    point3 = Array{SVector{2,Float64}}(undef, Nx, Ny*N_subs)
    point4 = Array{SVector{2,Float64}}(undef, Nx, Ny*N_subs)

    # Loop over all sub-enclosures.
    for k = 1:N_subs

        xs_now = xs[k] # Reference point.
        ys_now = ys[k]

        # Get change of coordinates at outer dimensions.
        # x first.
        deltaXbot = xs_now[2]-xs_now[1]
        deltaXtop = xs_now[4]-xs_now[3]
        deltaXleft = xs_now[5]-xs_now[4]
        # then y.
        deltaYright = ys_now[3]-ys_now[2]

        # Loop to create sub-divisions (cells) of each sub-enclosure.
        for m = 1:Ny+1
            # Add a bit of y.
            dy = (m-1)*deltaYright/Ny
            
            # We move the reference point in each iteration.
            refmoveXleft = (m-1)*deltaXleft/Ny
            refmoveXright = deltaXbot - (m-1)*(deltaXbot+deltaXtop)/Ny

            # Start at the bottom left point in any enclosure.
            for n = 1:Nx+1
                # Change in x and y.
                # First term = bottom left point.
                # Second term = how much we move the reference.
                # Third term = how much we move right in each iteration.
                xPoints[n,m] = xs_now[1] - refmoveXleft + (n-1)*(refmoveXright)/Nx
                yPoints[n,m] = ys_now[1] + dy
            end
        end

        # Put the points into the matrices to work with ray tracing.
        for m = 1:Ny
            yCount = m+(k-1)*Ny
            for n = 1:Nx
                point1[n,yCount] = SVector(xPoints[n,m], yPoints[n,m])
                point2[n,yCount] = SVector(xPoints[n+1,m], yPoints[n,m])
                point3[n,yCount] = SVector(xPoints[n+1,m+1], yPoints[n+1,m+1])
                point4[n,yCount] = SVector(xPoints[n,m+1], yPoints[n,m+1])
            end
        end
    end
    
    # Option to view geometry.
    if displayGeometry
        plot!(aspect_ratio=1.0)
        for m = 1:Nx # Plot imaginary walls between points to show the domain.
            for n = 1:N_subs*Ny
                # Plot all cells if requested.
                pointOne = point1[m,n] # Set the bottom left point.
                pointTwo = point2[m,n] # Set the bottom right point.
                pointThree = point3[m,n] # Set the top right point.
                pointFour = point4[m,n] # Set the top left point.
                display(plot!(trapezoid(pointOne, pointTwo, pointThree, pointFour),
                                legend=false,color=:white))
            end
        end
    end

    # miscellaneous useful numbers related to geometry
    N_surfs = 2*Nx+2*Ny*N_subs
    N_vols = Nx*Ny*N_subs

    return point1, point2, point3, point4, N_surfs, N_vols

end