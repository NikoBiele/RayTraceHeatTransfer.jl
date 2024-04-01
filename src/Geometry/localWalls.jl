"""
    localWalls(point1::Vector{Float64},point2::Vector{Float64},
                point3::Vector{Float64},point4::Vector{Float64})

This function determines a point on each wall of the current cell
as well as an inward pointing normal vector of each wall.
This information fully describe the walls of each cell.
"""
function localWalls(point1::Vector{Float64},point2::Vector{Float64},
                    point3::Vector{Float64},point4::Vector{Float64})

    # Points on each wall    
    # for output readability
    wallPointBottom = point1
    wallPointRight = point2
    wallPointTop = point3
    wallPointLeft = point4

    # bottom wall, unit vectors in local coordinate system
    xVecLocal = SVector{2}(normalize(point2-point1))
    yVecLocal = SVector{2}([-xVecLocal[2], xVecLocal[1]])
    bottomWallNormal = yVecLocal

    # right wall, unit vectors in local coordinate system
    xVecLocal = SVector{2}(normalize(point3-point2))
    yVecLocal = SVector{2}([-xVecLocal[2], xVecLocal[1]])
    rightWallNormal = yVecLocal

    # top wall, unit vectors in local coordinate system
    xVecLocal = SVector{2}(normalize(point4-point3))
    yVecLocal = SVector{2}([-xVecLocal[2], xVecLocal[1]])
    topWallNormal = yVecLocal

    # left wall, unit vectors in local coordinate system
    xVecLocal = SVector{2}(normalize(point1-point4))
    yVecLocal = SVector{2}([-xVecLocal[2], xVecLocal[1]])
    leftWallNormal = yVecLocal

    return wallPointBottom, wallPointRight, wallPointTop, wallPointLeft, bottomWallNormal, rightWallNormal, topWallNormal, leftWallNormal
end