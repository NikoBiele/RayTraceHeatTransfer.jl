"""
    SubEnclosure

This struct is used for initial geometry building.
To Initialize it use:
SubEnclosure(point1::Vector{Float64},point2::Vector{Float64},
                point3::Vector{Float64},point4::Vector{Float64},
                solidWall1::Bool,solidWall2::Bool,
                solidWall3::Bool,solidWall4::Bool)
The first four inputs are the bounding points.
The next four inputs state whether the walls are solid or empty.
"""
struct SubEnclosure
    # this struct is used when generating the mesh in a modular way

    point1::SVector{2,Float64} # bottom left vertex
    point2::SVector{2,Float64} # bottom right vertex
    point3::SVector{2,Float64} # top right vertex
    point4::SVector{2,Float64} # top left vertex
    solidWall1::Bool # bottom wall
    solidWall2::Bool # right wall
    solidWall3::Bool # top wall
    solidWall4::Bool # left wall
    midpoint::SVector{2,Float64} # midpoint of sub-enclosure
    wallPointBottom::SVector{2,Float64} # a point on the bottom wall
    wallPointRight::SVector{2,Float64} # a point on the right wall
    wallPointTop::SVector{2,Float64} # a point on the top wall
    wallPointLeft::SVector{2,Float64} # a point on the left wall
    bottomWallNormal::SVector{2,Float64} # inward facing normal of bottom wall
    rightWallNormal::SVector{2,Float64} # inward facing normal of right wall
    topWallNormal::SVector{2,Float64} # inward facing normal of top wall
    leftWallNormal::SVector{2,Float64} # inward facing normal of left wall
    function SubEnclosure(point1,point2,point3,point4,
                    solidWall1,solidWall2,solidWall3,solidWall4)
        wallPointBottom, wallPointRight, wallPointTop, wallPointLeft,
        bottomWallNormal, rightWallNormal, topWallNormal, leftWallNormal =
                    localWalls(point1,point2,point3,point4)
        return new(SVector(point1[1],point1[2]),
                    SVector(point2[1],point2[2]),
                    SVector(point3[1],point3[2]),
                    SVector(point4[1],point4[2]),
                    solidWall1,solidWall2,solidWall3,solidWall4,
                    SVector{2}((point1.+point2.+point3.+point4)./4),
                    wallPointBottom, wallPointRight, wallPointTop, wallPointLeft,
                    bottomWallNormal, rightWallNormal, topWallNormal, leftWallNormal)
    end
end