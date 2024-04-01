"""
    coefs_exchange(mesh::RayTracingMesh, Wall_absorbX::Array{Float64, 3},
                        Wall_absorbY::Array{Float64, 3},
                        N_abs_gas::Array{Float64, 3},RayCountTotal::Float64)

This function is called internallly during the ray tracing.
After ray tracing one emitter (surface or volume) this function calculates
one row of the two corresponding exchange factor matrices.
"""
function coefs_exchange(mesh::RayTracingMesh, Wall_absorbX::Array{Float64, 3}, Wall_absorbY::Array{Float64, 3},
                        N_abs_gas::Array{Float64, 3},RayCountTotal::Float64)
    
    # setup one row of FSS or FGS matrix, same function for both (rays absorbed by walls)
    FXS = zeros(1, mesh.N_surfs)
    wallCount = 0
    for subs = 1:mesh.N_subs
        for wall = 1:4
            if mesh.solidWalls[subs][wall]
                if wall == 1 || wall == 3 # bottom or top wall
                    if wall == 1
                        xRange = 1:mesh.Nx
                        yRange = 1
                    elseif wall == 3
                        xRange = 1:mesh.Nx
                        yRange = mesh.Ny
                    end
                    for m in xRange # loop through all surfaces
                        for n in yRange
                            wallCount +=1
                            FXS[1,wallCount] = Wall_absorbY[subs,m,n]/RayCountTotal
                        end
                    end
                else #if wall == 2 || wall == 4, left or right wall
                    if wall == 2
                        xRange = mesh.Nx
                        yRange = 1:mesh.Ny
                    elseif wall == 4
                        xRange = 1
                        yRange = 1:mesh.Ny
                    end
                    for m in xRange # loop through all surfaces
                        for n in yRange
                            wallCount +=1
                            FXS[1,wallCount] = Wall_absorbX[subs,m,n]/RayCountTotal
                        end
                    end
                end
            end
        end
    end

    # setup one row of FSG or FGG matrix, same function for both (rays absorbed by volumes)
    FXG = zeros(1, mesh.N_vols)
    volCount = 0
    for subs = 1:mesh.N_subs
        for m = 1:mesh.Nx
            for n = 1:mesh.Ny
                volCount += 1
                FXG[1,volCount] = N_abs_gas[subs,m,n]/RayCountTotal
            end
        end
    end

    return FXS, FXG
end