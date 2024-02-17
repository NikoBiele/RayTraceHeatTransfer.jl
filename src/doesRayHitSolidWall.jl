function doesRayHitSolidWall(point::SVector{2,Float64}, S::Float64, u_real::Float64, u_index::Int64, dir::SVector{2,Float64},
                            alphaWalls::Float64, Nx_coarse::Int64, Ny_coarse::Int64, xCount_coarse::Int64, yCount_coarse::Int64,
                            N_subs::Int64, point1_fine::Matrix{SVector{2,Float64}}, point2_fine::Matrix{SVector{2,Float64}},
                            point3_fine::Matrix{SVector{2,Float64}}, point4_fine::Matrix{SVector{2,Float64}},
                            Nx_fine::Int64, Ny_fine::Int64, logicalCores::Int64, NeighborIndices::Vector{SVector{2,Int64}},
                            Wall_absorbX::Array{Int64}, Wall_absorbY::Array{Int64},displayWhileTracing::Bool,
                            point1_coarse::Matrix{SVector{2,Float64}}, point2_coarse::Matrix{SVector{2,Float64}},
                            point3_coarse::Matrix{SVector{2,Float64}}, point4_coarse::Matrix{SVector{2,Float64}})

    # find out if any of the boundaries of the current cell are not allowed to be crossed (solid walls)
    solidWallx, solidWally = solidWall(Nx_coarse, Ny_coarse, xCount_coarse, yCount_coarse, N_subs)

    # then we found out if our ray hit a solid wall
    if u_index == solidWallx || u_index == solidWally # ray hits a solid wall

        pointOld = point # save old point before update (for plotting)
        S = S - u_real # update the travelled distance
        point = point + (u_real-1e-9)*dir # go almost into the wall (for plotting)

        # then go to reflect or absorp
        R_alpha = rand() # sample surface
        if R_alpha > alphaWalls # currently alpha == 1 so always absorbed (due to exchange factor approach)
            # then the ray is reflected diffusely
            # sample a new direction on the given surface
            # use a rotation matrix to transfer the angle from local to global

        else # then the ray is absorbed

            if displayWhileTracing
                display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = "",))
                display(scatter!((point[1], point[2]), color = "green", label = "", markersize = 2))
            end

            # Here we get the position in the fine grid
            xCount_fine, yCount_fine = whichEnclosure(point, point1_fine, point2_fine, point3_fine, point4_fine, N_subs, Nx_fine, Ny_fine)

            # increase a counter for the wall which was hit in order to measure the factors
            if u_index == solidWallx # necessary to discern between x and y for enclosures with 2 walls
                Wall_absorbX[xCount_fine,yCount_fine,logicalCores] += 1 # increase counter for the given wall
            elseif u_index == solidWally
                Wall_absorbY[xCount_fine,yCount_fine,logicalCores] += 1 # increase counter for the given wall
            else
                println("Unknown wall error in absorption.")
            end

            # then emit a new ray
            rayWasAbsorbed = true
            willRayHitWall = false
        end
    else # ray hits an imaginary boundary

        pointOld = point # save old point before update
        point = point + (u_real+1e-10).*dir # transfer ray to next enclosure
        
        if displayWhileTracing
            display(plot!([pointOld[1], point[1]], [pointOld[2], point[2]], label = ""))
        end

        S -= u_real # update the remaining distance

        u_real, u_index,xCount_coarse, yCount_coarse = distToNearestSurfAlongDir(point, dir, point1_coarse, point2_coarse, point3_coarse, point4_coarse,
                                                                    N_subs, Nx_coarse, Ny_coarse, NeighborIndices)

        # check if ray hits wall
        willRayHitWall = true
        rayWasAbsorbed = false
    end
    
    return rayWasAbsorbed, willRayHitWall, u_real, u_index, point, S, Wall_absorbX, Wall_absorbY, xCount_coarse, yCount_coarse
end