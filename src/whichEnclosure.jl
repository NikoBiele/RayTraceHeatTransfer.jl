function whichEnclosure(point::SVector{2,Float64}, mesh::TracingMesh, fineMesh::Bool)
    # this function finds out which enclosure we are in

    # rewrite the current location point
    xPoint = point[1]
    yPoint = point[2]
    xp = xPoint[1]
    yp = yPoint[1]

    # initiate outputs
    xCount = 0
    yCount = 0
    N_subs_count = 0.0

    # loop over all control volumes to find the one we are in
    # if neighbors is given as an input we only need to check those neighbors
    # else we check all cells
    if fineMesh
        for m = 1:mesh.Nx
            for n = 1:mesh.Ny*mesh.N_subs 
                Dcount = 0 # reset counters
                # get the four bounding points of this control volume
                @inbounds pointOne = mesh.point1[m,n] 
                @inbounds pointTwo = mesh.point2[m,n] 
                @inbounds pointThree = mesh.point3[m,n] 
                @inbounds pointFour = mesh.point4[m,n] 
                for p = 1:4 
                    if p == 1 
                        x1 = pointOne[1] 
                        y1 = pointOne[2] 
                        x2 = pointTwo[1] 
                        y2 = pointTwo[2] 
                    elseif p == 2 
                        x1 = pointTwo[1] 
                        y1 = pointTwo[2] 
                        x2 = pointThree[1] 
                        y2 = pointThree[2] 
                    elseif p == 3 
                        x1 = pointThree[1] 
                        y1 = pointThree[2] 
                        x2 = pointFour[1] 
                        y2 = pointFour[2] 
                    elseif p == 4 
                        x1 = pointFour[1] 
                        y1 = pointFour[2] 
                        x2 = pointOne[1] 
                        y2 = pointOne[2] 
                    end 
                    # check if the point is inside this control volume 
                    # check that it is to the left of each boundary 
                    D = (x2[1] - x1[1]) * (yp[1] - y1[1]) - (xp[1] - x1[1]) * (y2[1] - y1[1]) 
                    if D > 0 # then the point is on the right hand side 
                        Dcount += 1 
                    else # if point is on left hand side increment counter 
                        break   
                    end 
                    if Dcount == 4 # if point is to the left of all points 
                        xCount = m 
                        yCount = n
                        return xCount, yCount
                    end
                end
            end
        end
    else # if coarse 
        for j = 1:mesh.N_subs
            Dcount = 0 # reset counters
            # get the four bounding points of this control volume
            pointOne = mesh.point1_coarse[1,j] 
            pointTwo = mesh.point2_coarse[1,j] 
            pointThree = mesh.point3_coarse[1,j] 
            pointFour = mesh.point4_coarse[1,j] 
            for p = 1:4 
                if p == 1 
                    x1 = pointOne[1] 
                    y1 = pointOne[2] 
                    x2 = pointTwo[1] 
                    y2 = pointTwo[2] 
                elseif p == 2 
                    x1 = pointTwo[1] 
                    y1 = pointTwo[2] 
                    x2 = pointThree[1] 
                    y2 = pointThree[2] 
                elseif p == 3 
                    x1 = pointThree[1] 
                    y1 = pointThree[2] 
                    x2 = pointFour[1] 
                    y2 = pointFour[2] 
                elseif p == 4 
                    x1 = pointFour[1] 
                    y1 = pointFour[2] 
                    x2 = pointOne[1] 
                    y2 = pointOne[2] 
                end 
                # check if the point is inside this control volume 
                # check that it is to the left of each boundary 
                D = (x2[1] - x1[1]) * (yp[1] - y1[1]) - (xp[1] - x1[1]) * (y2[1] - y1[1]) 
                if D > 0 # then the point is on the right hand side 
                    Dcount += 1 
                else # if point is on left hand side increment counter 
                    break   
                end 
                if Dcount == 4 # if point is to the left of all points 
                    N_subs_count = j
                    return N_subs_count
                end
            end
        end
    end
end