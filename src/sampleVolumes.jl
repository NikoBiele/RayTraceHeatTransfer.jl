function sampleVolumes(point1_coarse::Matrix{SVector{2,Float64}}, point2_coarse::Matrix{SVector{2,Float64}},
                        point3_coarse::Matrix{SVector{2,Float64}}, point4_coarse::Matrix{SVector{2,Float64}}, Ny_coarse::Int, Nx_coarse::Int,
                        N_surfs_fine::Int,N_vols_fine::Int,point1_fine::Matrix{SVector{2,Float64}}, point2_fine::Matrix{SVector{2,Float64}},
                        point3_fine::Matrix{SVector{2,Float64}}, point4_fine::Matrix{SVector{2,Float64}}, Ny_fine::Int, Nx_fine::Int,
                        beta::Float64,omega::Float64,N_rays::Int64,displayWhileTracing::Bool,nthreads::Int,N_subs::Int,
                        NeighborIndices_coarse::Matrix{Array{SVector{2, Int64}}})
    
    # This function loops over all of the volumetric (gas) emitters.
    # Each sampled emitter generates rows in the FGS and FGG matrices.

    println("Starting volume sampling:")

    @time begin
        # Initialize the matrices.
        FGS = zeros(N_vols_fine, N_surfs_fine)
        FGG = zeros(N_vols_fine, N_vols_fine)

        # Tell the program that we are not sampling surfaces.
        wallEmitter = 0
        sampleLeftRight = false
        sampleTopBottom = false

        # Calculate number of gas volume emitters to be sampled.
        N_vols = Nx_fine*Ny_fine*N_subs
        volumeEmitter = 0

        # Loop over all volumetric emitters.
        for i = 1:Nx_fine
            for j = 1:Ny_fine*N_subs
                volumeEmitter += 1
                println("Now sampling volume emitter number $volumeEmitter/$N_vols.") # progress update in the REPL.
                xCountSample = i
                yCountSample = j

                # Call ray tracing function and sample current volume.
                Wall_absorbX, Wall_absorbY, N_abs_gas, RayCountTotal = rayTracing_2D(point1_coarse, point2_coarse,
                                                                                    point3_coarse, point4_coarse,
                                                                                    Ny_coarse, Nx_coarse, point1_fine,
                                                                                    point2_fine, point3_fine, point4_fine,
                                                                                    Ny_fine, Nx_fine, beta,
                                                                                    omega, N_rays, displayWhileTracing,
                                                                                    nthreads, wallEmitter, volumeEmitter,
                                                                                    N_subs, xCountSample, yCountSample,
                                                                                    sampleLeftRight, sampleTopBottom,
                                                                                    NeighborIndices_coarse)
                # Determine coefficients in matrices.
                FGS[volumeEmitter, :], FGG[volumeEmitter, :] = coefs_FGS_FGG(Wall_absorbX, Wall_absorbY, N_abs_gas, N_surfs_fine, 
                                                                        N_vols_fine, RayCountTotal, Nx_fine, Ny_fine, N_subs)    
            end
        end
        println("Volume sampling finished, the elapsed time for volume sampling was:")
    end

    return FGS, FGG
    
end