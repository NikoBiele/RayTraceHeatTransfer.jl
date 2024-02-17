function doesRayAbsorbOrScatter(point::SVector{2,Float64}, pointOld::SVector{2,Float64}, point1_fine::Matrix{SVector{2,Float64}},
                                point2_fine::Matrix{SVector{2,Float64}}, point3_fine::Matrix{SVector{2,Float64}},
                                point4_fine::Matrix{SVector{2,Float64}}, N_subs::Int64, Nx_fine::Int64, Ny_fine::Int64,
                                logicalCores::Int64, omega::Float64, N_abs_gas::Array{Int64})
    # initialize outputs
    dir = SVector(0.0, 0.0)
    S = 0.0

    # if we absorp, we increase counter and emit a new ray
    # if we scatter we sample a new direction from a scattering phase function
    R_omega = rand()
    if R_omega > omega
        # ray is absorbed, increase counter and emit a new ray

        # here we find out which (fine) enclosure we are in
        xCount_fine, yCount_fine = whichEnclosure(point, point1_fine, point2_fine, point3_fine, point4_fine, N_subs, Nx_fine, Ny_fine)

        N_abs_gas[xCount_fine, yCount_fine, logicalCores] += 1 # increase counter (on fine grid)

        if displayWhileTracing # green for absorption
            display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = ""))
            display(scatter!((point[1], point[2]), color = "green", label = "", markersize = 2))
        end

        rayWasAbsorbed = true

    else # R_omega < omega
        # ray is scattered

        if displayWhileTracing # red for scattering
            display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = ""))
            display(scatter!((point[1], point[2]), color = "red", label = "", markersize = 2))
        end

        # Find isotropic scatter direction (3D spherical projected onto 2D)
        theta = acos(2.0*rand()-1.0) # cone angle
        phi = 2*pi*rand() # circumferential angle (turn the cone angle around the plane)
        # direction vector according to conventions
        r = 1
        xdir = r*sin(theta)*cos(phi)
        # ydir = r*sin(theta)*sin(phi) # this one is projected onto plane (set to zero)
        zdir = r*cos(theta) 
        dir = SVector(xdir, zdir)

        # now we emit a new ray from the scatter point
        R_S = rand() # sample for finding distance travelled by ray
        S = -(1/beta)*log(R_S) # distance travelled by ray

        rayWasAbsorbed = false
    end

    return rayWasAbsorbed, dir, S, N_abs_gas
end