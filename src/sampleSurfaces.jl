function sampleSurfaces(mesh::TracingMesh,gas::GasProperties,N_rays::Int64,nthreads::Int64,displayWhileTracing::Bool)
    
    # This function loops over all of the surface emitters.
    # Each sampled emitter generates rows in the FSS and FSG matrices.

    println("Starting surface sampling:")

    @time begin
        # Initialize the matrices.
        FSS = zeros(mesh.N_surfs, mesh.N_surfs)
        FSG = zeros(mesh.N_surfs, mesh.N_vols)

        # Tell the program we are sampling surfaces, not volumes
        volumeEmitter = 0 # set to zero for no volume emission
        sampleLeftRight = true # sample left-right walls of enclosure
        sampleTopBottom = false # do not sample top-bottom of enclosure
        wallEmitter = 0 # surface counter
        N_walls = 2*mesh.Nx+2*mesh.Ny*mesh.N_subs

        for m = [1, mesh.Nx] # sample left wall, then right wall
            for n = 1:mesh.Ny*mesh.N_subs # from bottom to top
                wallEmitter += 1
                println("Now sampling wall emitter number $wallEmitter/$N_walls.") # progress update in the REPL
                
                # set x- and y-counters
                xCountSample = m
                yCountSample = n

                # sample left-right walls one by one
                Wall_absorbX_sum, Wall_absorbY_sum, N_abs_gas_sum, RayCountTotal_sum = rayTracing_2D(
                                                                                mesh, gas, N_rays,
                                                                                displayWhileTracing, nthreads,
                                                                                wallEmitter, volumeEmitter,
                                                                                xCountSample, yCountSample,
                                                                                sampleLeftRight, sampleTopBottom)
                # determine coefficients in matrices (to save the ray tracing)
                FSS[wallEmitter, :], FSG[wallEmitter, :] = coefs_FSS_FSG(mesh, Wall_absorbX_sum, Wall_absorbY_sum,
                                                                            N_abs_gas_sum, RayCountTotal_sum)  
            end
        end

        # sample top and bottom
        sampleLeftRight = false
        sampleTopBottom = true

        for n = [1, mesh.Ny*mesh.N_subs] # sample bottom, then top
            for m = 1:mesh.Nx # sample along the width
                wallEmitter += 1
                println("Now sampling wall emitter number $wallEmitter/$N_walls.") # progress update in the REPL

                # set x- and y-counters
                xCountSample = m
                yCountSample = n

                # sample top-bottom walls one by one
                Wall_absorbX_sum, Wall_absorbY_sum, N_abs_gas_sum, RayCountTotal_sum = rayTracing_2D( 
                                                                        mesh, gas, N_rays,
                                                                        displayWhileTracing, nthreads,
                                                                        wallEmitter, volumeEmitter,
                                                                        xCountSample, yCountSample,
                                                                        sampleLeftRight, sampleTopBottom)
                # determine coefficients in matrices (to save the ray tracing)
                FSS[wallEmitter, :], FSG[wallEmitter, :] = coefs_FSS_FSG(mesh, Wall_absorbX_sum, Wall_absorbY_sum,
                                                                            N_abs_gas_sum, RayCountTotal_sum)    
            end   
        end
        println("Surface sampling finished, the elapsed time for surface sampling was:")
    end

    return FSS, FSG
    
end