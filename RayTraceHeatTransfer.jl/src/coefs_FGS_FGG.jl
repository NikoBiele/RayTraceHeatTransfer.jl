function coefs_FGS_FGG(Wall_absorbX::Matrix{Float64}, Wall_absorbY::Matrix{Float64},
                        N_abs_gas::Matrix{Float64}, N_surfs, N_vols, RayCountTotal::Float64, Nx::Int, Ny::Int, N_subs::Int)
    # This function calculates one row of the FGS matrix and one row of the FGG-matrix
    # this calculation is based on the ray tracing results of one volume emitter

    FGS = zeros(1, N_surfs)
    for m = [1, Nx] # only take ends
        for n = 1:Ny*N_subs # loop through all positions
            if m == 1 # left hand
                FGS[1,n] = Wall_absorbX[m,n]/RayCountTotal
            elseif m == Nx # right side
                FGS[1,n+Ny*N_subs] = Wall_absorbX[m,n]/RayCountTotal
            end
        end
    end
    for m = 1:Nx # loop through all positions
        for n = [1, Ny*N_subs] # only need top and bottom
            if n == 1 # bottom
                FGS[1,2*Ny*N_subs+m] = Wall_absorbY[m,n]/RayCountTotal
            elseif n == Ny*N_subs # top
                FGS[1,2*Ny*N_subs+Nx+m] = Wall_absorbY[m,n]/RayCountTotal
            end
        end
    end
    # setup one row of FGG matrix
    FGG = zeros(1, N_vols)
    volCount = 0
    for m = 1:Nx
        for n = 1:Ny*N_subs
            volCount += 1
            FGG[1,volCount] = N_abs_gas[m,n]/RayCountTotal
        end
    end

    return FGS, FGG
end