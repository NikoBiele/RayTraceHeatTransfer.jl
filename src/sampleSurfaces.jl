function sampleSurfaces(point1_coarse::Matrix{SVector{2,Float64}}, point2_coarse::Matrix{SVector{2,Float64}},
                        point3_coarse::Matrix{SVector{2,Float64}}, point4_coarse::Matrix{SVector{2,Float64}}, Ny_coarse::Int, Nx_coarse::Int,
                        N_surfs_fine::Int,N_vols_fine::Int,point1_fine::Matrix{SVector{2,Float64}}, point2_fine::Matrix{SVector{2,Float64}},
                        point3_fine::Matrix{SVector{2,Float64}}, point4_fine::Matrix{SVector{2,Float64}}, Ny_fine::Int, Nx_fine::Int,
                        beta::Float64,omega::Float64,N_rays::Int64,displayWhileTracing::Bool,nthreads::Int,N_subs::Int,
                        NeighborIndices_coarse::Matrix{Array{SVector{2, Int64}}})
    
    # This function loops over all of the surface emitters.
    # Each sampled emitter generates rows in the FSS and FSG matrices.

    println("Starting surface sampling:")

    @time begin
        # Initialize the matrices.
        FSS = zeros(N_surfs_fine, N_surfs_fine)
        FSG = zeros(N_surfs_fine, N_vols_fine)

        # Tell the program we are sampling surfaces, not volumes
        volumeEmitter = 0 # set to zero for no volume emission
        sampleLeftRight = true # sample left-right walls of enclosure
        sampleTopBottom = false # do not sample top-bottom of enclosure
        wallEmitter = 0 # surface counter
        N_walls = 2*Nx_fine+2*Ny_fine

        for m = [1, Nx_fine] # sample left wall, then right wall
            for n = 1:Ny_fine*N_subs # from bottom to top
                wallEmitter += 1
                println("Now sampling wall emitter number $wallEmitter/$N_walls.") # progress update in the REPL
                
                # set x- and y-counters
                xCountSample = m
                yCountSample = n

                # sample left-right walls one by one
                Wall_absorbX, Wall_absorbY, N_abs_gas, RayCountTotal = rayTracing_2D(
                                                            point1_coarse, point2_coarse, point3_coarse, 
                                                            point4_coarse, Ny_coarse, Nx_coarse,
                                                            point1_fine, point2_fine, point3_fine,
                                                            point4_fine, Ny_fine, Nx_fine,
                                                            beta, omega, N_rays,
                                                            displayWhileTracing, nthreads,
                                                            wallEmitter, volumeEmitter, N_subs,
                                                            xCountSample, yCountSample,
                                                            sampleLeftRight, sampleTopBottom,
                                                            NeighborIndices_coarse)

                # determine coefficients in matrices (to save the ray tracing)
                FSS[wallEmitter, :], FSG[wallEmitter, :] = coefs_FSS_FSG(Wall_absorbX, Wall_absorbY,
                                                                        N_abs_gas, N_surfs_fine, N_vols_fine,
                                                                        RayCountTotal, Nx_fine, Ny_fine, N_subs)  
            end
        end

        # sample top and bottom
        sampleLeftRight = false
        sampleTopBottom = true

        for n = [1, Ny_fine*N_subs] # sample bottom, then top
            for m = 1:Nx_fine # sample along the width
                wallEmitter += 1
                println("Now sampling wall emitter number $wallEmitter/$N_walls.") # progress update in the REPL

                # set x- and y-counters
                xCountSample = m
                yCountSample = n

                # sample top-bottom walls one by one
                Wall_absorbX, Wall_absorbY, N_abs_gas, RayCountTotal = rayTracing_2D( 
                                                                    point1_coarse, point2_coarse, point3_coarse, 
                                                                    point4_coarse, Ny_coarse, Nx_coarse,
                                                                    point1_fine, point2_fine, point3_fine,
                                                                    point4_fine, Ny_fine, Nx_fine,
                                                                    beta, omega, N_rays,
                                                                    displayWhileTracing, nthreads,
                                                                    wallEmitter, volumeEmitter, N_subs,
                                                                    xCountSample, yCountSample,
                                                                    sampleLeftRight, sampleTopBottom,
                                                                    NeighborIndices_coarse)
                # determine coefficients in matrices (to save the ray tracing)
                FSS[wallEmitter, :], FSG[wallEmitter, :] = coefs_FSS_FSG(Wall_absorbX, Wall_absorbY,
                                                                    N_abs_gas, N_surfs_fine, N_vols_fine,
                                                                    RayCountTotal, Nx_fine, Ny_fine, N_subs)    
            end   
        end
        println("Surface sampling finished, the elapsed time for surface sampling was:")
    end

    return FSS, FSG
    
end