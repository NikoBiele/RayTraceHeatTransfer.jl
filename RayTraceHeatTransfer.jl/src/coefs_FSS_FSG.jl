function coefs_FSS_FSG(Wall_absorbX::Matrix{Float64}, Wall_absorbY::Matrix{Float64},
                        N_abs_gas::Matrix{Float64}, N_surfs::Int64, N_vols::Int64,
                        RayCountTotal::Float64, Nx::Int, Ny::Int, N_subs::Int)
    # This function calculates one row of the FSS matrix and one row of the FSG-matrix
    # this calculation is based on the ray tracing results of one wall emitter
    
    # setup one row of FSS matrix
    FSS = zeros(1, N_surfs)
    for m = [1, Nx] # only take ends
        for n = 1:Ny*N_subs # loop through all surfaces
            if m == 1 # left hand
                FSS[1,n] = Wall_absorbX[m,n]/RayCountTotal
            elseif m == Nx # right side
                FSS[1,n+Ny*N_subs] = Wall_absorbX[m,n]/RayCountTotal
            end
        end
    end
    for m = 1:Nx # loop through all positions
        for n = [1, Ny*N_subs] # only need top and bottom
            if n == 1 # bottom
                FSS[1,2*Ny*N_subs+m] = Wall_absorbY[m,n]/RayCountTotal
            elseif n == Ny*N_subs # top
                FSS[1,2*Ny*N_subs+Nx+m] = Wall_absorbY[m,n]/RayCountTotal
            end
        end
    end

    # setup one row of FSG matrix
    FSG = zeros(1, N_vols)
    volCount = 0
    for m = 1:Nx
        for n = 1:Ny*N_subs
            volCount += 1
            FSG[1,volCount] = N_abs_gas[m,n]/RayCountTotal
        end
    end
    return FSS, FSG
end