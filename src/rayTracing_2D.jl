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
    N_abs_gas = Array{Int64}(undef, Nx_fine, Ny_fine*N_subs, nthreads) # zeros(Nx_fine, Ny_fine*N_subs, nthreads) # set gas absorption counters to zero
    Wall_absorbX = Array{Int64}(undef, Nx_fine, Ny_fine*N_subs, nthreads) # zeros(Nx_fine, Ny_fine*N_subs, nthreads) # set wall absorption counters to zero
    Wall_absorbY = Array{Int64}(undef, Nx_fine, Ny_fine*N_subs, nthreads) # zeros(Nx_fine, Ny_fine*N_subs, nthreads) # set wall absorption counters to zero
    RayCountTotal = Vector{Int64}(undef, nthreads) # zeros(nthreads)  # counters for how many rays have been emitted from emissive power
    N_abs_gas .= 0
    Wall_absorbX .= 0
    Wall_absorbY .= 0    
    RayCountTotal .= 0
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
                point, dir = sampleSurface(Nx_fine, Ny_fine, xCountSample, yCountSample, N_subs,
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
                point, dir = sampleVolume(Nx_fine, Ny_fine, N_subs, volumeEmitter,
                                            point1_fine, point2_fine, point3_fine, point4_fine,xCountSample,yCountSample)

            else
                println("Unknown volume emitter.")
            end
            
            # from here the emission point and direction is known

            if displayWhileTracing # plot emission point (yellow for emission)
                display(scatter!((point[1], point[2]), color = "yellow", label = "", markersize = 2))
            end

            R_S = rand() # sample for finding ray distance
            S = -(1/beta)*log(R_S) # get the ray distance travelled before absorption or scattering

            # now we ray trace the first ray in the enclosure in which the emission happens
            # first we find out where in the coarse mesh we are
            # as well as the distance to the nearest wall

            u_real, u_index, xCount_coarse, yCount_coarse = distToNearestSurfAlongDir(point, dir, point1_coarse, point2_coarse, point3_coarse, point4_coarse,
                                                                                        N_subs, Nx_coarse, Ny_coarse)

            # now we determine if the ray hit a wall
            @label doesRayHitWall
            if S < u_real
                # S < u_real, so ray does not hit a wall but is scattered or absorbed in gas
                # bundle is scattered or absorbed in the participating medium
                while S < u_real
                    # this loop keeps running as long as the ray keeps scattering
                    # this loop ends if the ray is absorbed in the gas or it hits a wall

                    # update position and save old point for plotting
                    pointOld = point
                    point = point + (S-1e-9)*dir

                    # means we hit the gas and either absorp or scatters
                    rayWasAbsorbed, dir, S = doesRayAbsorbOrScatter(point, pointOld, point1_fine, point2_fine, point3_fine, point4_fine,
                                                                    N_subs, Nx_fine, Ny_fine, logicalCores,omega,N_abs_gas,displayWhileTracing)
                    if rayWasAbsorbed
                        @goto emitNewRay
                    else
                        u_real, u_index, xCount_coarse, yCount_coarse = distToNearestSurfAlongDir(point, dir, point1_coarse, point2_coarse, point3_coarse, point4_coarse,
                                                                                N_subs, Nx_coarse, Ny_coarse, NeighborIndices_coarse[xCount_coarse,yCount_coarse])
                    end
                end
                # we go out of the loop above when the remaining distance is longer than the distance to the wall
                # therefore we break out of the loop above when we hit a wall!

                # find out if the wall we hit was solid
                rayWasAbsorbed, willRayHitWall, u_real, u_index, point, S, Wall_absorbX, Wall_absorbY, xCount_coarse, yCount_coarse = 
                                                                doesRayHitSolidWall(point, S, u_real, u_index, dir, alphaWalls, Nx_coarse, Ny_coarse,
                                                                                    xCount_coarse, yCount_coarse, N_subs, point1_fine, point2_fine,
                                                                                    point3_fine, point4_fine, Nx_fine, Ny_fine, logicalCores,
                                                                                    NeighborIndices_coarse[xCount_coarse,yCount_coarse],
                                                                                    Wall_absorbX, Wall_absorbY,displayWhileTracing,
                                                                                    point1_coarse, point2_coarse, point3_coarse, point4_coarse)
                if rayWasAbsorbed
                    @goto emitNewRay
                elseif willRayHitWall
                    @goto doesRayHitWall
                end
            else
                # u_real < S, so bundle hits a wall! (real or imaginary)

                # find out if the wall we hit was solid
                rayWasAbsorbed, willRayHitWall, u_real, u_index, point, S, Wall_absorbX, Wall_absorbY, xCount_coarse, yCount_coarse = 
                                                                doesRayHitSolidWall(point, S, u_real, u_index, dir, alphaWalls, Nx_coarse, Ny_coarse,
                                                                                    xCount_coarse, yCount_coarse, N_subs, point1_fine, point2_fine,
                                                                                    point3_fine, point4_fine, Nx_fine, Ny_fine, logicalCores,
                                                                                    NeighborIndices_coarse[xCount_coarse,yCount_coarse],
                                                                                    Wall_absorbX, Wall_absorbY,displayWhileTracing,
                                                                                    point1_coarse, point2_coarse, point3_coarse, point4_coarse)
                if rayWasAbsorbed
                    @goto emitNewRay
                elseif willRayHitWall
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