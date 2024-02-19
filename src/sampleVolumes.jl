function sampleVolumes(mesh::TracingMesh,gas::GasProperties,N_rays::Int64,nthreads::Int64,displayWhileTracing::Bool)
    
    # This function loops over all of the volumetric (gas) emitters.
    # Each sampled emitter generates rows in the FGS and FGG matrices.

    println("Starting volume sampling:")

    @time begin
        # Initialize the matrices.
        FGS = zeros(mesh.N_vols, mesh.N_surfs)
        FGG = zeros(mesh.N_vols, mesh.N_vols)

        # Tell the program that we are not sampling surfaces.
        wallEmitter = 0
        sampleLeftRight = false
        sampleTopBottom = false

        # Calculate number of gas volume emitters to be sampled.
        N_vols = mesh.Nx*mesh.Ny*mesh.N_subs
        volumeEmitter = 0

        # Loop over all volumetric emitters.
        for i = 1:mesh.Nx
            for j = 1:mesh.Ny*mesh.N_subs
                volumeEmitter += 1
                println("Now sampling volume emitter number $volumeEmitter/$N_vols.") # progress update in the REPL.

                xCountSample = i
                yCountSample = j

                # Call ray tracing function and sample current volume.
                Wall_absorbX_sum, Wall_absorbY_sum, N_abs_gas_sum, RayCountTotal_sum = rayTracing_2D(mesh, gas, N_rays,
                                                                                    displayWhileTracing, nthreads,
                                                                                    wallEmitter, volumeEmitter, xCountSample, yCountSample,
                                                                                    sampleLeftRight, sampleTopBottom)
                # Determine coefficients in matrices.
                FGS[volumeEmitter, :], FGG[volumeEmitter, :] = coefs_FGS_FGG(mesh, Wall_absorbX_sum, Wall_absorbY_sum,
                                                                                N_abs_gas_sum, RayCountTotal_sum)    
            end
        end
        println("Volume sampling finished, the elapsed time for volume sampling was:")
    end

    return FGS, FGG
    
end