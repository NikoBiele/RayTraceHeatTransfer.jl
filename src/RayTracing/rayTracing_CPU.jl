"""
    rayTracing_CPU(subNumber::Int64,wallNumber::Int64,
                        sampleLeftRight::Bool,sampleTopBottom::Bool,
                        mesh::RayTracingMesh, gas::GasProperties, N_rays::Int64,
                        wallEmission::Bool,volumeEmission::Bool,
                        displayWhileTracing::Bool,nthreads::Int64,
                        xCountSample::Int64, yCountSample::Int64)

This function ray traces a single ray in a 2D domain (or a number of rays in parallel).
The ray (or rays) is emitted from a specified emitter (wall or volume) and traced throughout the domain.
This function also gives the user the possibility to plot the traced ray paths (must execute single threaded).
When all the rays have been traced the absorbed rays are summed up for each cell
and returned for use in the calculation of the exchange factor matrices.
The ray tracing is performed on a coarse grid and absorption points are mapped to a fine grid.
This approach greatly improves the efficiency as opposed to tracing on a fine grid.
"""
function rayTracing_CPU(subNumber::Int64,wallNumber::Int64,
                        sampleLeftRight::Bool,sampleTopBottom::Bool,
                        mesh::RayTracingMesh, gas::GasProperties, N_rays::Int64,
                        wallEmission::Bool,volumeEmission::Bool,
                        displayWhileTracing::Bool,nthreads::Int64,
                        xCountSample::Int64, yCountSample::Int64)
  
    # set counters to zero for each logical core
    N_abs_gas = zeros(mesh.N_subs, mesh.Nx, mesh.Ny, nthreads) # set gas absorption counters to zero
    Wall_absorbX = zeros(mesh.N_subs, mesh.Nx, mesh.Ny, nthreads) # set wall absorption counters to zero
    Wall_absorbY = zeros(mesh.N_subs, mesh.Nx, mesh.Ny, nthreads) # set wall absorption counters to zero
    RayCountTotal = zeros(nthreads)  # counters for how many rays have been emitted from emissive power
    rays_per_thread = trunc(Int, round(N_rays/nthreads)) # total number of rays to emit

    # set absorptivity of walls (walls always absorb)
    alphaWalls = 1.0 # must be = 1

    # enter the multi-threaded loop
    Threads.@threads for logicalCores = 1:nthreads
 
        # first check if all of the rays have been emitted
        while RayCountTotal[logicalCores] <= rays_per_thread

            point = [] # reset the emission point
            RayCountTotal[logicalCores] += 1 # increment emission counter

            # emitter is specified as a surface
            if wallEmission
                # sample surface
                @inline point, i1 = sampleSurface(mesh, subNumber, wallNumber, xCountSample, yCountSample, sampleLeftRight, sampleTopBottom)
            end

            # emitter is specified as a volume
            if volumeEmission
                # sample volume
                @inline point, i1 = sampleVolume(mesh, xCountSample, yCountSample)
            end
            
            # from here the emission point and direction is known

            if displayWhileTracing # plot emission point (yellow for emission)
                display(scatter!((point[1], point[2]), color = "yellow", label = "", markersize = 2))
            end

            # here we find out which (coarse) enclosure we are in
            @inline N_subs_count = whichSubEnclosure(point, mesh)

            R_S = rand() # sample for finding ray distance
            S = -(1/gas.beta)*log(R_S) # get the ray distance travelled before absorption or scattering
            
            # now we ray trace the first ray in the enclosure in which the emission happens    
            # we calculate the point on each wall based on the control volume we are currently in
            # get distance to the surface along this direction
            @inline u_real, u_index = distToSurface(point, i1, mesh, N_subs_count)

            # now we determine if the ray hit a wall
            emitNewRay = false # ensure that we don't break out unintended
            doesRayHitWall = true
            while doesRayHitWall == true
                if S < u_real
                    # S < u_real, so ray does not hit a wall but is scattered or absorbed in gas
                    # bundle is scattered or absorbed in the participating medium
                    while S < u_real
                        # this loop keeps running as long as the ray keeps scattering
                        # this loop ends if the ray is absorbed in the gas or it hits a wall

                        # update position and save old point for plotting
                        pointOld = point
                        point = point + (S-1e-12)*i1

                        # means we hit the gas and either absorp or scatters
                        # if we absorp, we increase counter and emit a new ray
                        # if we scatter we sample a new direction from a scattering phase function
                        R_omega = rand()
                        if R_omega > gas.omega
                            # ray is absorbed, increase counter and emit a new ray

                            # here we find out which (fine) enclosure we are in
                            @inline N_subs_count = whichSubEnclosure(point, mesh)
                            @inline xCount, yCount = whichCell(point, mesh)
                            yCount = yCount - mesh.Ny*(N_subs_count-1)

                            N_abs_gas[N_subs_count, xCount, yCount, logicalCores] += 1 # increase counter (on fine grid)

                            if displayWhileTracing # green for absorption
                                display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = ""))
                                display(scatter!((point[1], point[2]), color = "green", label = "", markersize = 2))
                            end

                            emitNewRay = true
                            break

                        else # R_omega < omega
                            # ray is scattered

                            if displayWhileTracing # red for scattering
                                display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = ""))
                                display(scatter!((point[1], point[2]), color = "red", label = "", markersize = 2))
                            end

                            # sample an isotropic scattering direction
                            i1 = isotropicScatter()

                            # now we emit a new ray from the scatter point
                            R_S = rand() # sample for finding distance travelled by ray
                            S = -(1/gas.beta)*log(R_S) # distance travelled by ray
                            
                            # here we find out which (coarse) enclosure we are
                            # fineMesh = false
                            @inline N_subs_count = whichSubEnclosure(point, mesh)

                            # get distances to (coarse) surfaces
                            @inline u_real, u_index = distToSurface(point, i1, mesh, N_subs_count) 
                        end

                    end

                    # if emitNewRay was set to true, also break out of the doesRayHitWall-loop
                    if emitNewRay == true
                        break
                    end

                    # we go out of the loop above when the remaining distance is longer than the distance to the wall
                    # therefore we break out of the loop above when we hit a wall!

                    # find out if any of the boundaries of the current cell are not allowed to be crossed (solid walls)
                    @inline solidWall_1, solidWall_2, solidWall_3, solidWall_4 = solidWalls(mesh, N_subs_count)             
                    wallCount_1 = 0
                    wallCount_2 = 0
                    wallCount_3 = 0
                    wallCount_4 = 0
                    if solidWall_1
                        wallCount_1 = 1
                    end
                    if solidWall_2
                        wallCount_2 = 2
                    end
                    if solidWall_3
                        wallCount_3 = 3
                    end       
                    if solidWall_4
                        wallCount_4 = 4
                    end               

                    # then we found out if our ray hit a solid wall
                    if u_index == wallCount_1 || u_index == wallCount_2 || u_index == wallCount_3 || u_index == wallCount_4 # ray hits a solid wall

                        pointOld = point # save old point before update (for plotting)
                        S = S - u_real # update the travelled distance
                        point = point + (u_real-1e-12)*i1 # go almost into the wall (for plotting)

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
                            # fineMesh = true
                            # @inline N_subs_count = whichSubEnclosure(point, mesh)
                            @inline xCount, yCount = whichCell(point, mesh) #, N_subs_count)

                            # increase a counter for the wall which was hit in order to measure the factors
                            if u_index == wallCount_2 || u_index == wallCount_4 # necessary to discern between x and y for enclosures with 2 walls
                                Wall_absorbX[N_subs_count,xCount,yCount,logicalCores] += 1 # increase counter for the given wall
                            elseif u_index == wallCount_1 || u_index == wallCount_3
                                Wall_absorbY[N_subs_count,xCount,yCount,logicalCores] += 1 # increase counter for the given wall
                            else
                                println("Unknown wall error in absorption.")
                            end

                            # then emit a new ray
                            emitNewRay = true
                            break

                        end
                    else # ray hits an imaginary boundary

                        pointOld = point # save old point before update
                        point = point + (u_real+1e-12).*i1 # transfer ray to next enclosure
                        
                        if displayWhileTracing
                            display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = ""))
                        end

                        S -= u_real # update the remaining distance

                        # here we find out which (coarse) enclosure we ended up in
                        @inline N_subs_count = whichSubEnclosure(point, mesh)

                        if N_subs_count === nothing
                            # ray has escaped the domain (should be very rare)
                            # go back one small step and absorb the ray
                            # then break to emit new ray
                            point = point - 2*1e-12.*i1 # go back one small step (into nearby cell)

                            # here we find out which (fine) enclosure we are in
                            @inline N_subs_count = whichSubEnclosure(point, mesh)
                            @inline xCount, yCount = whichCell(point, mesh)
                            yCount = yCount - mesh.Ny*(N_subs_count-1)

                            N_abs_gas[N_subs_count, xCount, yCount, logicalCores] += 1 # increase counter (on fine grid)

                            if displayWhileTracing # green for absorption
                                display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = ""))
                                display(scatter!((point[1], point[2]), color = "green", label = "", markersize = 2))
                            end

                            emitNewRay = true
                            break
                        end

                        # get distance to surfaces
                        @inline u_real, u_index = distToSurface(point, i1, mesh, N_subs_count)
                
                        # check if ray hits wall
                        doesRayHitWall = true

                    end

                else
                    # u_real < S, so bundle hits a wall! (real or imaginary)

                    # find out if any of the boundaries of the current cell are not allowed to be crossed (solid walls)
                    @inline solidWall_1, solidWall_2, solidWall_3, solidWall_4 = solidWalls(mesh, N_subs_count)
                    wallCount_1 = 0
                    wallCount_2 = 0
                    wallCount_3 = 0
                    wallCount_4 = 0
                    if solidWall_1
                        wallCount_1 = 1
                    end
                    if solidWall_2
                        wallCount_2 = 2
                    end
                    if solidWall_3
                        wallCount_3 = 3
                    end       
                    if solidWall_4
                        wallCount_4 = 4
                    end

                    # then we found out if our ray hit a solid wall
                    if u_index == wallCount_1 || u_index == wallCount_2 || u_index == wallCount_3 || u_index == wallCount_4 # ray hits a solid wall
                        
                        pointOld = point # save old point before update (for plotting)
                        S = S - u_real # update the travelled distance
                        point = point + (u_real-1e-12)*i1 # go almost into the wall (for plotting)

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
                            # fineMesh = true
                            # @inline N_subs_count = whichSubEnclosure(point, mesh)
                            @inline xCount, yCount = whichCell(point, mesh) #, N_subs_count)
                            yCount = yCount - mesh.Ny*(N_subs_count-1)

                            # increase a counter for the wall which was hit in order to measure the factors
                            if u_index == wallCount_2 || u_index == wallCount_4 # necessary to discern between x and y for enclosures with 2 walls
                                Wall_absorbX[N_subs_count,xCount,yCount,logicalCores] += 1 # increase counter for the given wall
                            elseif u_index == wallCount_1 || u_index == wallCount_3
                                Wall_absorbY[N_subs_count,xCount,yCount,logicalCores] += 1 # increase counter for the given wall
                            else
                                println("Unknown wall error in absorption.")
                            end

                            # then emit a new ray
                            emitNewRay = true
                            break

                        end
                    else # ray hits an imaginary boundary

                        pointOld = point # save old point before update
                        point = point + (u_real+1e-12).*i1 # transfer ray to next enclosure

                        if displayWhileTracing
                            display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = ""))
                        end

                        S -= u_real # update the remaining distance

                        # here we find out which (coarse) enclosure we ended up in
                        # fineMesh = false
                        @inline N_subs_count = whichSubEnclosure(point, mesh)

                        if N_subs_count === nothing
                            # ray has escaped the domain (should be very rare)
                            # go back one small step and absorb the ray
                            # then break to emit new ray
                            point = point - 2*1e-12.*i1 # go back one small step (into nearby cell)

                            # here we find out which (fine) enclosure we are in
                            @inline N_subs_count = whichSubEnclosure(point, mesh)
                            @inline xCount, yCount = whichCell(point, mesh)
                            yCount = yCount - mesh.Ny*(N_subs_count-1)

                            N_abs_gas[N_subs_count, xCount, yCount, logicalCores] += 1 # increase counter (on fine grid)

                            if displayWhileTracing # green for absorption
                                display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = ""))
                                display(scatter!((point[1], point[2]), color = "green", label = "", markersize = 2))
                            end

                            emitNewRay = true
                            break
                        end
                        
                        # get distance to surfaces
                        @inline u_real, u_index = distToSurface(point, i1, mesh, N_subs_count)  
                                                        
                        # check if ray hits wall
                        doesRayHitWall = true
                    end
                end
            end
        end
    end

    # sum up the outputs over all threads
    Wall_absorbX_sum = dropdims(sum(Wall_absorbX, dims=4),dims=4)
    Wall_absorbY_sum = dropdims(sum(Wall_absorbY, dims=4),dims=4)
    N_abs_gas_sum = dropdims(sum(N_abs_gas, dims=4),dims=4)
    RayCountTotal_sum = sum(RayCountTotal)

    return Wall_absorbX_sum, Wall_absorbY_sum, N_abs_gas_sum, RayCountTotal_sum
end