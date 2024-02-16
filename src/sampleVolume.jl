function sampleVolume(Nx, Ny, N_subs, volumeEmitter::Int, point1::Matrix{SVector{2,Float64}},
                        point2::Matrix{SVector{2,Float64}}, point3::Matrix{SVector{2,Float64}}, point4::Matrix{SVector{2,Float64}},
                        xCount::Int64,yCount::Int64)
    
    # This function samples an emission point and direction in a volume (a surface in 2D), by separating it in
    # two triangles and sampling according to their area.
    # This operation requires knowledge of the mesh, which is why the point matrices are also inputs.
    
    # get the vertices of the cell
    A = point1[xCount, yCount]
    B = point2[xCount, yCount]
    C = point3[xCount, yCount]
    D = point4[xCount, yCount]

    # split into two triangles
    # calculate the area of the triangles
    ABC_area = 0.5*(A[1]*(B[2]-C[2])+B[1]*(C[2]-A[2])+C[1]*(A[2]-B[2]))
    CDA_area = 0.5*(C[1]*(D[2]-A[2])+D[1]*(A[2]-C[2])+A[1]*(C[2]-D[2]))

    # sample to determine which triangle to choose
    if rand() < ABC_area/(ABC_area+CDA_area)
        # sample ABC triangle
        R_1 = rand()
        R_2 = rand()
        sqrt_R1 = sqrt(R_1)
        point = (1-sqrt_R1)*A + sqrt_R1*(1-R_2)*B + sqrt_R1*R_2*C
    else
        # sample CDA triangle
        R_1 = rand()
        R_2 = rand()
        sqrt_R1 = sqrt(R_1)
        point = (1-sqrt_R1)*C + sqrt_R1*(1-R_2)*D + sqrt_R1*R_2*A
    end

    # now get direction
    # Sample a 3D spherical sample and project it onto the 2D plane (omit third component)
    theta = acos(2.0*rand()-1.0) # cone angle
    phi = 2*pi*rand() # circumferential angle (turn the cone angle around the plane)
    # direction vector according to conventions
    r = 1
    xdir = r*sin(theta)*cos(phi)
    # ydir = r*sin(theta)*sin(phi) # this one is projected onto plane (set to zero) 
    zdir = r*cos(theta)
    i1 = SVector(xdir, zdir)

    return point, i1
end