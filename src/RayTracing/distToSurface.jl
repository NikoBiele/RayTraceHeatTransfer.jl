"""
    distToSurface(point::SVector{2,Float64}, i1::SVector{2,Float64},
                    mesh::RayTracingMesh, N_subs_count::Int64)    

This function determines the distance to a plane (a line in 2D) from a
given point along a given direction, as well as the index of the wall (bottom, right, top, left).
"""
function distToSurface(point::SVector{2,Float64}, i1::SVector{2,Float64},
                        mesh::RayTracingMesh, N_subs_count::Int64)                   

    # lookup information about the local walls
    wallPointBottom = mesh.wallPointBottom[N_subs_count]
    wallPointRight = mesh.wallPointRight[N_subs_count]
    wallPointTop = mesh.wallPointTop[N_subs_count]
    wallPointLeft = mesh.wallPointLeft[N_subs_count]
    bottomWallNormal = mesh.bottomWallNormal[N_subs_count]
    rightWallNormal = mesh.rightWallNormal[N_subs_count]
    topWallNormal = mesh.topWallNormal[N_subs_count]
    leftWallNormal = mesh.leftWallNormal[N_subs_count]

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