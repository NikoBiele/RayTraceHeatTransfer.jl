"""
    sampleSurfaces(mesh::RayTracingMesh,gas::GasProperties,
                    N_rays::Int64,nthreads::Int64,displayWhileTracing::Bool)

This function samples all of the surface emitters.
Each sampled emitter generates rows in the FSS and FSG exchange factor matrices.                
"""
function sampleSurfaces(mesh::RayTracingMesh,gas::GasProperties,
                    N_rays::Int64,nthreads::Int64,displayWhileTracing::Bool)

    println("Starting surface sampling:")

    @time begin
        # Initialize the matrices.
        FSS = zeros(mesh.N_surfs, mesh.N_surfs)
        FSG = zeros(mesh.N_surfs, mesh.N_vols)
        # Tell the program we are sampling surfaces, not volumes
        wallEmission = true
        volumeEmission = false # set to zero for no volume emission
        N_walls = sum(sum(mesh.solidWalls))*mesh.Nx
        xCountSample = 0 # these are for volume sampling
        yCountSample = 0 # these are for volume sampling
        sampleLeftRight = false
        sampleTopBottom = false
        wallEmitter = 0 # surface counter

        for subCount = 1:mesh.N_subs

            # first sample bottom of each sub-enclosure (if it has a wall)
            sampleLeftRight = false
            sampleTopBottom = true
            j = 1
            wallCount = 1
    
            if mesh.solidWalls[subCount][wallCount] # if this wall is solid, sample it
                for i = 1:mesh.Nx

                    wallEmitter += 1
                    println("Now sampling wall emitter number $wallEmitter/$N_walls.") # progress update in the REPL

                    xCountSample = i
                    yCountSample = j

                    # sample left-right walls one by one
                    Wall_absorbX_sum, Wall_absorbY_sum, N_abs_gas_sum, RayCountTotal_sum = rayTracing_CPU(subCount,wallCount,
                                                                                    sampleLeftRight,sampleTopBottom,
                                                                                    mesh, gas, N_rays, wallEmission, volumeEmission,
                                                                                    displayWhileTracing, nthreads, xCountSample, yCountSample)
                    # determine coefficients in matrices (to save the ray tracing)
                    FSS[wallEmitter, :], FSG[wallEmitter, :] = coefs_exchange(mesh, Wall_absorbX_sum, Wall_absorbY_sum,
                                                                                N_abs_gas_sum, RayCountTotal_sum)  
                end
            end

            # next, sample right side
            sampleLeftRight = true
            sampleTopBottom = false
            i = mesh.Nx
            wallCount = 2
    
            if mesh.solidWalls[subCount][wallCount] # if this wall is solid, sample it
                for j = 1:mesh.Ny

                    wallEmitter += 1
                    println("Now sampling wall emitter number $wallEmitter/$N_walls.") # progress update in the REPL

                    xCountSample = i
                    yCountSample = j

                    # sample left-right walls one by one
                    Wall_absorbX_sum, Wall_absorbY_sum, N_abs_gas_sum, RayCountTotal_sum = rayTracing_CPU(subCount,wallCount,
                                                                                    sampleLeftRight,sampleTopBottom,
                                                                                    mesh, gas, N_rays, wallEmission, volumeEmission,
                                                                                    displayWhileTracing, nthreads, xCountSample, yCountSample)
                    # determine coefficients in matrices (to save the ray tracing)
                    FSS[wallEmitter, :], FSG[wallEmitter, :] = coefs_exchange(mesh, Wall_absorbX_sum, Wall_absorbY_sum,
                                                                                N_abs_gas_sum, RayCountTotal_sum)  
                end
            end
            
            # sample the top wall
            sampleLeftRight = false
            sampleTopBottom = true
            j = mesh.Ny
            wallCount = 3
        
            if mesh.solidWalls[subCount][wallCount] # if this wall is solid, sample it
                for i = 1:mesh.Nx

                    wallEmitter += 1
                    println("Now sampling wall emitter number $wallEmitter/$N_walls.") # progress update in the REPL

                    xCountSample = i
                    yCountSample = j

                    # sample left-right walls one by one
                    Wall_absorbX_sum, Wall_absorbY_sum, N_abs_gas_sum, RayCountTotal_sum = rayTracing_CPU(subCount,wallCount,
                                                                                    sampleLeftRight,sampleTopBottom,
                                                                                    mesh, gas, N_rays, wallEmission, volumeEmission,
                                                                                    displayWhileTracing, nthreads, xCountSample, yCountSample)
                    # determine coefficients in matrices (to save the ray tracing)
                    FSS[wallEmitter, :], FSG[wallEmitter, :] = coefs_exchange(mesh, Wall_absorbX_sum, Wall_absorbY_sum,
                                                                                N_abs_gas_sum, RayCountTotal_sum)  
                end
            end

            sampleLeftRight = true
            sampleTopBottom = false
            i = 1
            wallCount = 4
    
            if mesh.solidWalls[subCount][wallCount] # if this wall is solid, sample it
                for j = 1:mesh.Ny
    
                    wallEmitter += 1
                    println("Now sampling wall emitter number $wallEmitter/$N_walls.") # progress update in the REPL

                    xCountSample = i
                    yCountSample = j

                    # sample left-right walls one by one
                    Wall_absorbX_sum, Wall_absorbY_sum, N_abs_gas_sum, RayCountTotal_sum = rayTracing_CPU(subCount,wallCount,
                                                                                    sampleLeftRight,sampleTopBottom,
                                                                                    mesh, gas, N_rays, wallEmission, volumeEmission,
                                                                                    displayWhileTracing, nthreads, xCountSample, yCountSample)
                    # determine coefficients in matrices (to save the ray tracing)
                    FSS[wallEmitter, :], FSG[wallEmitter, :] = coefs_exchange(mesh, Wall_absorbX_sum, Wall_absorbY_sum,
                                                                                N_abs_gas_sum, RayCountTotal_sum)  
                end
            end
        end
        println("Surface sampling finished, the elapsed time for surface sampling was:")
    end

    return FSS, FSG
    
end