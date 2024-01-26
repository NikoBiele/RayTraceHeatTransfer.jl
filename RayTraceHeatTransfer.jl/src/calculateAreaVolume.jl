function calculateAreaVolume(Nx::Int64,Ny::Int64,N_subs::Int64,
                            width::Float64,point1,point2,point3,point4)

    # this function calculates the area of walls and the volume of volumes
    # First calculate the area of walls
    # calculate the length of lines and multiply by a width

    A_count = 0
    Area = zeros(2*Nx+2*Ny*N_subs)
    for m = 1 # left wall
        for n = 1:Ny*N_subs
            pointOne = point1[m,n] # extract point
            pointFour = point4[m,n]
            dx = abs(pointOne[1]-pointFour[1])
            dy = abs(pointOne[2]-pointFour[2])
            A_count += 1 # increment counter
            Area[A_count] = sqrt(dx^2+dy^2)*width # calculate area
        end
    end
    for m = Nx # right wall
        for n = 1:Ny*N_subs
            pointTwo = point2[m,n] # extract point
            pointThree = point3[m,n]
            dx = abs(pointTwo[1]-pointThree[1])
            dy = abs(pointTwo[2]-pointThree[2])
            A_count += 1 # increment counter
            Area[A_count] = sqrt(dx^2+dy^2)*width # calculate area
        end
    end
    for m = 1:Nx # bottom wall
        for n = 1
            pointOne = point1[m,n] # extract point
            pointTwo = point2[m,n]
            dx = abs(pointOne[1]-pointTwo[1])
            dy = abs(pointOne[2]-pointTwo[2])
            A_count += 1 # increment counter
            Area[A_count] = sqrt(dx^2+dy^2)*width # calculate area
        end
    end
    for m = 1:Nx # top wall
        for n = Ny*N_subs
            pointThree = point3[m,n] # extract point
            pointFour = point4[m,n]
            dx = abs(pointThree[1]-pointFour[1])
            dy = abs(pointThree[2]-pointFour[2])
            A_count += 1 # increment counter
            Area[A_count] = sqrt(dx^2+dy^2)*width # calculate area
        end
    end
    # calculate the volume of volumes
    V_count = 0
    Volume = zeros(Nx*Ny*N_subs)
    for m = 1:Nx
        for n = 1:Ny*N_subs
            # get the vertices of the enclosure
            A = point1[m, n]
            B = point2[m, n]
            C = point3[m, n]
            D = point4[m, n]
            # split into two triangles
            # calculate the area of the triangles
            ABC_area = 0.5*(A[1]*(B[2]-C[2])+B[1]*(C[2]-A[2])+C[1]*(A[2]-B[2]))
            CDA_area = 0.5*(C[1]*(D[2]-A[2])+D[1]*(A[2]-C[2])+A[1]*(C[2]-D[2]))
            V_count += 1
            Area_tot = ABC_area + CDA_area
            Volume[V_count] = Area_tot*width
        end
    end
    
    return Area, Volume
end