function distToNearestSurfAlongDir(point::SVector{2,Float64}, dir::SVector{2,Float64}, point1_coarse::Matrix{SVector{2,Float64}},
                                    point2_coarse::Matrix{SVector{2,Float64}}, point3_coarse::Matrix{SVector{2,Float64}},
                                    point4_coarse::Matrix{SVector{2,Float64}}, N_subs::Int64, Nx_coarse::Int64, Ny_coarse::Int64, NeighborIndices=nothing)

    # here we find out which (coarse) enclosure we are
    xCount_coarse, yCount_coarse = whichEnclosure(point, point1_coarse, point2_coarse, point3_coarse, point4_coarse,
                                                    N_subs, Nx_coarse, Ny_coarse, NeighborIndices)

    # we calculate the point on each wall based on the control volume we are currently in
    wallPointBottom, wallPointRight, wallPointTop, wallPointLeft, bottomWallNormal, rightWallNormal, topWallNormal, leftWallNormal =
                                    localWalls(xCount_coarse, yCount_coarse, point1_coarse, point2_coarse, point3_coarse, point4_coarse)

    # get distances to (coarse) surfaces
    u_real, u_index = distToSurface(point, dir, wallPointTop, wallPointBottom, wallPointLeft, wallPointRight,
                                    topWallNormal, bottomWallNormal, leftWallNormal, rightWallNormal)

    return u_real, u_index, xCount_coarse, yCount_coarse

end