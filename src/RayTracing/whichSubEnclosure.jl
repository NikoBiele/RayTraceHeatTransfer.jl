"""
    whichSubEnclosure(point::SVector{2,Float64},mesh::RayTracingMesh)

This function finds out which SubEnclosure the ray is currently in.
"""
function whichSubEnclosure(point::SVector{2,Float64},mesh::RayTracingMesh)

    # rewrite the current location point
    xPoint = point[1]
    yPoint = point[2]
    xp = xPoint[1]
    yp = yPoint[1]

    # initiate counter
    N_subs_count = 0

    for j = 1:mesh.N_subs
        Dcount = 0 # reset counters
        # get the four bounding points of this control volume
        @inbounds pointOne = mesh.point1_coarse[1,j]
        @inbounds pointTwo = mesh.point2_coarse[1,j]
        @inbounds pointThree = mesh.point3_coarse[1,j]
        @inbounds pointFour = mesh.point4_coarse[1,j]
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

    return nothing
end