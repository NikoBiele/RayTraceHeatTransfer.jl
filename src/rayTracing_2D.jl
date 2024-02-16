function rayTracing_2D(point1_coarse::Matrix{SVector{2,Float64}}, point2_coarse::Matrix{SVector{2,Float64}}, point3_coarse::Matrix{SVector{2,Float64}}, 
                        point4_coarse::Matrix{SVector{2,Float64}}, Ny_coarse::Int, Nx_coarse::Int,
                        point1_fine::Matrix{SVector{2,Float64}}, point2_fine::Matrix{SVector{2,Float64}}, point3_fine::Matrix{SVector{2,Float64}},
                        point4_fine::Matrix{SVector{2,Float64}}, Ny_fine::Int, Nx_fine::Int,
                        beta::Float64, omega::Float64, N_rays::Int64,
                        displayWhileTracing::Bool, nthreads::Int64,
                        wallEmitter::Int, volumeEmitter::Int, N_subs::Int, xCountSample::Int, yCountSample::Int,
                        sampleLeftRight::Bool, sampleTopBottom::Bool,
                        NeighborIndices_coarse::Matrix{Array{SVector{2, Int64}}})

    # This function ray traces a single ray in a 2D domain (or a number of rays in parallel).
    # The ray (or rays) is emitted from a specified emitter (wall or volume) and
    # traced throughout the domain.
    # This function also gives the user the possibility to plot the traced ray paths (must execute single threaded).
    # When all the rays have been traced the absorbed rays are summed up for each cell
    # and returned for use in the exchange factor matrices.
    # The ray tracing is performed on a coarse grid and absorption points are mapped to a fine grid.
    # This approach greatly improves the efficiency as opposed to tracing on a fine grid.

    # set counters to zero for each logical core
    N_abs_gas = zeros(Nx_fine, Ny_fine*N_subs, nthreads) # set gas absorption counters to zero
    Wall_absorbX = zeros(Nx_fine, Ny_fine*N_subs, nthreads) # set wall absorption counters to zero
    Wall_absorbY = zeros(Nx_fine, Ny_fine*N_subs, nthreads) # set wall absorption counters to zero
    RayCountTotal = zeros(nthreads)  # counters for how many rays have been emitted from emissive power
    rays_per_thread = trunc(Int, round(N_rays/nthreads)) # total number of rays to emit
    
    # set absorptivity of walls (walls always absorb)
    alphaWalls = 1.0 # must be = 1

    # enter the multi-threaded loop
    Threads.@threads for logicalCores = 1:nthreads
 
        @label emitNewRay # emit a new ray
        point = [] # reset the emission point
        RayCountTotal[logicalCores] += 1 # increment emission counter

        # first check if all of the rays have been emitted
        if RayCountTotal[logicalCores] <= rays_per_thread

            # emitter is specified as a surface
            if wallEmitter == 0 # set to zero if we do not emit from any surface
                # empty on purpose, set to zero for no emission
                # to be explicit about when we are sampling a volume

            elseif wallEmitter >= 1 && wallEmitter <= 2*Nx_fine+2*Ny_fine*N_subs

                # sample surface
                point, i1 = sampleSurface(Nx_fine, Ny_fine, xCountSample, yCountSample, N_subs,
                                            sampleLeftRight, sampleTopBottom, point1_fine, point2_fine,
                                            point3_fine, point4_fine)

            else
                println("Unknown surface emitter.")
            end

            # emitter is specified as a volume
            if volumeEmitter == 0 # set to zero if we do not emit from any volume
                # left empty on purpose, set to zero for no emission
                # to be explicit about when we are sampling a surface

            elseif volumeEmitter >= 1 && volumeEmitter <= Nx_fine*Ny_fine*N_subs 

                # sample volume
                point, i1 = sampleVolume(Nx_fine, Ny_fine, N_subs, volumeEmitter,
                                            point1_fine, point2_fine, point3_fine, point4_fine,xCountSample,yCountSample)

            else
                println("Unknown volume emitter.")
            end
            
            # from here the emission point and direction is known

            if displayWhileTracing # plot emission point (yellow for emission)
                display(scatter!((point[1], point[2]), color = "yellow", label = "", markersize = 2))
            end

            # here we find out which (coarse) enclosure we are in
            xCount_coarse, yCount_coarse = whichEnclosure(point, point1_coarse, point2_coarse, point3_coarse, point4_coarse,
                                                               N_subs, Nx_coarse, Ny_coarse)

            R_S = rand() # sample for finding ray distance
            S = -(1/beta)*log(R_S) # get the ray distance travelled before absorption or scattering
            
            # now we ray trace the first ray in the enclosure in which the emission happens
            # we do this by calculating the distance to all walls of this enclosure
            # first we define the walls
            # for this we need a point on each wall and a normal vector
            
            # we calculate the point on each wall based on the control volume we are currently in
            wallPointBottom, wallPointRight, wallPointTop, wallPointLeft, bottomWallNormal, rightWallNormal, topWallNormal, leftWallNormal =
                                            localWalls(xCount_coarse, yCount_coarse, point1_coarse, point2_coarse, point3_coarse, point4_coarse)

            # get distance to the surface along this direction
            u_real, u_index = distToSurface(point, i1, wallPointTop, wallPointBottom, wallPointLeft, wallPointRight,
                                                    topWallNormal, bottomWallNormal, leftWallNormal, rightWallNormal)

            # now we determine if the ray hit a wall
            @label doesRayHitWall
            if S < u_real
                # S < u_real, so ray does not hit a wall but is scattered or absorbed in gas
                @label scatter_or_absorb # bundle is scattered or absorbed in the participating medium

                # update position and save old point for plotting
                pointOld = point
                point = point + (S-1e-9)*i1

                # means we hit the gas and either absorp or scatters
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

                    @goto emitNewRay

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
                    i1 = SVector(xdir, zdir)

                end

                # now we emit a new ray from the scatter point
                R_S = rand() # sample for finding distance travelled by ray
                S = -(1/beta)*log(R_S) # distance travelled by ray
                
                # here we find out which (coarse) enclosure we are
                xCount_coarse, yCount_coarse = whichEnclosure(point, point1_coarse, point2_coarse, point3_coarse, point4_coarse,
                                                                    N_subs, Nx_coarse, Ny_coarse, NeighborIndices_coarse[xCount_coarse,yCount_coarse])

                # we calculate the point on each wall based on the control volume we are currently in
                wallPointBottom, wallPointRight, wallPointTop, wallPointLeft, bottomWallNormal, rightWallNormal, topWallNormal, leftWallNormal =
                                                    localWalls(xCount_coarse, yCount_coarse, point1_coarse, point2_coarse, point3_coarse, point4_coarse)

                # get distances to (coarse) surfaces
                u_real, u_index = distToSurface(point, i1, wallPointTop, wallPointBottom, wallPointLeft, wallPointRight, topWallNormal, bottomWallNormal, leftWallNormal, rightWallNormal)

                if u_real > S
                    # ray is scattered or absorbed
                    @goto scatter_or_absorb

                else
                    # u_real < S, so bundle hits a wall! (real or imaginary)
                    @goto rayHitsWall

                end
            else
                @label rayHitsWall
                # u_real < S, so bundle hits a wall! (real or imaginary)

                # find out if any of the boundaries of the current cell are not allowed to be crossed (solid walls)
                solidWallx, solidWally = solidWall(Nx_coarse, Ny_coarse, xCount_coarse, yCount_coarse, N_subs)

                # then we found out if our ray hit a solid wall
                if u_index == solidWallx || u_index == solidWally # ray hits a solid wall

                    pointOld = point # save old point before update (for plotting)
                    S = S - u_real # update the travelled distance
                    point = point + (u_real-1e-9)*i1 # go almost into the wall (for plotting)

                    # then go to reflect or absorp
                    R_alpha = rand() # sample surface
                    if R_alpha > alphaWalls # currently alpha == 1 so always absorbed (due to exchange factor approach)
                        # then the ray is reflected diffusely
                        # sample a new direction on the given surface
                        # use a rotation matrix to transfer the angle from local to global

                    else # then the ray is absorbed

                        if displayWhileTracing
                            display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = "",))
                            display(scatter!((point[1], point[2]), color = "green", label = "", markersize = 2))
                        end

                        # Here we get the position in the fine grid
                        xCount_fine, yCount_fine = whichEnclosure(point, point1_fine, point2_fine, point3_fine, point4_fine, N_subs, Nx_fine, Ny_fine)

                        # increase a counter for the wall which was hit in order to measure the factors
                        if u_index == solidWallx # necessary to discern between x and y for enclosures with 2 walls
                            Wall_absorbX[xCount_fine,yCount_fine,logicalCores] += 1 # increase counter for the given wall
                        elseif u_index == solidWally
                            Wall_absorbY[xCount_fine,yCount_fine,logicalCores] += 1 # increase counter for the given wall
                        else
                            println("Unknown wall error in absorption.")
                        end

                        # then emit a new ray
                        @goto emitNewRay

                    end
                else # ray hits an imaginary boundary

                    pointOld = point # save old point before update
                    point = point + (u_real+1e-10).*i1 # transfer ray to next enclosure
                    
                    if displayWhileTracing
                        display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = ""))
                    end

                    S -= u_real # update the remaining distance

                    # here we find out which (coarse) enclosure we ended up in
                    xCount_coarse, yCount_coarse = whichEnclosure(point, point1_coarse, point2_coarse, point3_coarse, point4_coarse,
                                                            N_subs, Nx_coarse, Ny_coarse,NeighborIndices_coarse[xCount_coarse,yCount_coarse])

                    # we calculate the point on each wall based on the control volume we are currently in
                    wallPointBottom, wallPointRight, wallPointTop, wallPointLeft, bottomWallNormal, rightWallNormal, topWallNormal, leftWallNormal =
                                                    localWalls(xCount_coarse, yCount_coarse, point1_coarse, point2_coarse, point3_coarse, point4_coarse)

                    # get distance to surfaces
                    u_real, u_index = distToSurface(point, i1, wallPointTop, wallPointBottom, wallPointLeft, 
                                                    wallPointRight, topWallNormal, bottomWallNormal, leftWallNormal, rightWallNormal)
            
                    # check if ray hits wall
                    @goto doesRayHitWall
                end
            end
        end
    end

    # sum up the outputs
    Wall_absorbX_tot = dropdims(sum(Wall_absorbX, dims=3),dims=3)
    Wall_absorbY_tot = dropdims(sum(Wall_absorbY, dims=3),dims=3)
    N_abs_gas_tot = dropdims(sum(N_abs_gas, dims=3),dims=3)
    RayCountTotal_tot = sum(RayCountTotal)

    return Wall_absorbX_tot, Wall_absorbY_tot, N_abs_gas_tot, RayCountTotal_tot
end