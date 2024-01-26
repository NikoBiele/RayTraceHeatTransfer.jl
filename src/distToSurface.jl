function distToSurface(point::SVector{2,Float64}, i1::SVector{2,Float64}, wallPointTop::SVector{2,Float64},
                        wallPointBottom::SVector{2,Float64}, wallPointLeft::SVector{2,Float64},
                        wallPointRight::SVector{2,Float64}, topWallNormal::SVector{2,Float64},
                        bottomWallNormal::SVector{2,Float64}, leftWallNormal::SVector{2,Float64}, rightWallNormal::SVector{2,Float64})
                        
    # This function determines the distance to a plane (a line in 2D) from a
    # given point along a given direction, as well as the index of the wall (bottom, right, top, left).

    # determine distance to all four walls from point along direction
    # some of these will be negative
    u_top = (dot((wallPointTop.-point),topWallNormal))/(dot(i1,topWallNormal))
    u_bottom = (dot((wallPointBottom.-point),bottomWallNormal))/(dot(i1,bottomWallNormal))
    u_left = (dot((wallPointLeft.-point),leftWallNormal))/(dot(i1,leftWallNormal))
    u_right = (dot((wallPointRight.-point),rightWallNormal))/(dot(i1,rightWallNormal))

    # find smallest positive distance and index of it, so we know which surface is hit
    u = [u_bottom, u_right, u_top, u_left]
    u_testForNegative = findall(x -> x <=0, u)
    u[u_testForNegative] .= Inf # set negative distances to infinity
    u_min_and_index = findmin(u) # then find minimum distance and its index
    u_real = u_min_and_index[1]
    u_index = u_min_and_index[2]

    return u_real, u_index
    
end