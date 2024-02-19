function coefs_FSS_FSG(mesh::TracingMesh, Wall_absorbX::Matrix{Float64}, Wall_absorbY::Matrix{Float64},
                        N_abs_gas::Matrix{Float64},RayCountTotal::Float64)
    # This function calculates one row of the FSS matrix and one row of the FSG-matrix
    # this calculation is based on the ray tracing results of one wall emitter
    
    # setup one row of FSS matrix
    FSS = zeros(1, mesh.N_surfs)
    for m = [1, mesh.Nx] # only take ends
        for n = 1:mesh.Ny*mesh.N_subs # loop through all surfaces
            if m == 1 # left hand
                FSS[1,n] = Wall_absorbX[m,n]/RayCountTotal
            elseif m == mesh.Nx # right side
                FSS[1,n+mesh.Ny*mesh.N_subs] = Wall_absorbX[m,n]/RayCountTotal
            end
        end
    end
    for m = 1:mesh.Nx # loop through all positions
        for n = [1, mesh.Ny*mesh.N_subs] # only need top and bottom
            if n == 1 # bottom
                FSS[1,2*mesh.Ny*mesh.N_subs+m] = Wall_absorbY[m,n]/RayCountTotal
            elseif n == mesh.Ny*mesh.N_subs # top
                FSS[1,2*mesh.Ny*mesh.N_subs+mesh.Nx+m] = Wall_absorbY[m,n]/RayCountTotal
            end
        end
    end

    # setup one row of FSG matrix
    FSG = zeros(1, mesh.N_vols)
    volCount = 0
    for m = 1:mesh.Nx
        for n = 1:mesh.Ny*mesh.N_subs
            volCount += 1
            FSG[1,volCount] = N_abs_gas[m,n]/RayCountTotal
        end
    end
    return FSS, FSG
end