"""
    areaVolumeMesh(Nx::Int64,Ny::Int64,N_subs::Int64,
                point1::Matrix{SVector{2,Float64}},point2::Matrix{SVector{2,Float64}},
                point3::Matrix{SVector{2,Float64}},point4::Matrix{SVector{2,Float64}})

This function is called internallly inside the RayTracingMesh struct when generating a new instance.
This function calculates the area of bounding walls and the volume of volumes.
"""
function areaVolumeMesh(Nx::Int64,Ny::Int64,N_subs::Int64,
                            point1::Matrix{SVector{2,Float64}},
                            point2::Matrix{SVector{2,Float64}},
                            point3::Matrix{SVector{2,Float64}},
                            point4::Matrix{SVector{2,Float64}})

    width = 1.0 # width of domain, third dimension (not important in 2D)
    Area = zeros(N_subs,2*Nx+2*Ny)
    Volume = zeros(N_subs,Nx,Ny)

    for i = 1:N_subs # for each sub enclosure
        wallCount = 0
        for m = 1:Nx # bottom wall
            for n = 1
                wallCount += 1 # increment counter
                # calculate area
                Area[i,wallCount] = norm(point1[m,n+(i-1)*Ny]-point2[m,n+(i-1)*Ny])*width
            end
        end
        for m = Nx # right wall
            for n = 1:Ny
                wallCount += 1 # increment counter
                # calculate area
                Area[i,wallCount] = norm(point2[m,n+(i-1)*Ny]-point3[m,n+(i-1)*Ny])*width
            end
        end
        for m = 1:Nx # top wall
            for n = Ny
                wallCount += 1 # increment counter
                # calculate area
                Area[i,wallCount] = norm(point3[m,n+(i-1)*Ny]-point4[m,n+(i-1)*Ny])*width
            end
        end
        for m = 1 # left wall
            for n = 1:Ny
                wallCount += 1 # increment counter
                # calculate area
                Area[i,wallCount] = norm(point1[m,n+(i-1)*Ny]-point4[m,n+(i-1)*Ny])*width
            end
        end

        # calculate the volume of volumes
        for m = 1:Nx
            for n = 1:Ny
                # get the vertices of the enclosure
                A = point1[m, n+(i-1)*Ny]
                B = point2[m, n+(i-1)*Ny]
                C = point3[m, n+(i-1)*Ny]
                D = point4[m, n+(i-1)*Ny]
                # split into two triangles
                # calculate the area of the triangles
                ABC_area = 0.5*(A[1]*(B[2]-C[2])+B[1]*(C[2]-A[2])+C[1]*(A[2]-B[2]))
                CDA_area = 0.5*(C[1]*(D[2]-A[2])+D[1]*(A[2]-C[2])+A[1]*(C[2]-D[2]))
                Area_tot = ABC_area + CDA_area
                Volume[i,m,n] = Area_tot*width
            end
        end
    end

    return Area, Volume
end