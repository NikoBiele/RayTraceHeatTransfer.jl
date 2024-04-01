"""
    RayTracingMesh

This struct generates the mesh from a vector of SubEnclosures when instantiated.
"""
struct RayTracingMesh

    # number of sub-enclosures
    N_subs::Int64
    subMidpoints::Vector{SVector{2,Float64}}

    # the coarse mesh
    Nx_coarse::Int64
    Ny_coarse::Int64
    point1_coarse::Matrix{SVector{2,Float64}}
    point2_coarse::Matrix{SVector{2,Float64}}
    point3_coarse::Matrix{SVector{2,Float64}}
    point4_coarse::Matrix{SVector{2,Float64}}
    N_surfs_coarse::Int64
    N_vols_coarse::Int64
    solidWalls::Vector{Vector{Bool}}
    # information about each coarse wall
    wallPointBottom::Vector{SVector{2,Float64}} # a point on the bottom wall
    wallPointRight::Vector{SVector{2,Float64}} # a point on the right wall
    wallPointTop::Vector{SVector{2,Float64}} # a point on the top wall
    wallPointLeft::Vector{SVector{2,Float64}} # a point on the left wall
    bottomWallNormal::Vector{SVector{2,Float64}} # inward facing normal of bottom wall
    rightWallNormal::Vector{SVector{2,Float64}} # inward facing normal of right wall
    topWallNormal::Vector{SVector{2,Float64}} # inward facing normal of top wall
    leftWallNormal::Vector{SVector{2,Float64}} # inward facing normal of left wall

    # the fine mesh
    Nx::Int64
    Ny::Int64
    point1::Matrix{SVector{2,Float64}}
    point2::Matrix{SVector{2,Float64}}
    point3::Matrix{SVector{2,Float64}}
    point4::Matrix{SVector{2,Float64}}
    N_surfs::Int64
    N_vols::Int64
    Area::Matrix{Float64}
    Volume::Array{Float64}

    # generate the fields automatically from fewer inputs about the bounding geometry
    function RayTracingMesh(subs::Vector{SubEnclosure},Ndim::Int64)
        # number of sub-enclosures
        N_subs = length(subs)
        # define the number of coarse splits in each enclosure
        # these are not allowed to be user defined
        Nx_coarse = 1
        Ny_coarse = 1
        # mesh to get coarse geometry
        point1_coarse, point2_coarse, point3_coarse, point4_coarse, N_surfs_coarse, N_vols_coarse =
                                                        meshGeometry(subs,Nx_coarse,Ny_coarse);
        # mesh to get fine geometry
        point1, point2, point3, point4, N_surfs, N_vols = meshGeometry(subs,Ndim,Ndim);
        # define the outer solid walls
        solidWalls = defineSolidWalls(subs)
        # calculate the area and volume of each cell in the fine mesh
        Area, Volume = areaVolumeMesh(Ndim,Ndim,N_subs,point1,point2,point3,point4);
        # find midpoint of each sub enclosure
        midPoints = [subs[i].midpoint for i = 1:N_subs]
        # setup wall information
        wallPointBottom = [subs[i].wallPointBottom for i = 1:N_subs]
        wallPointRight = [subs[i].wallPointRight for i = 1:N_subs]
        wallPointTop = [subs[i].wallPointTop for i = 1:N_subs]
        wallPointLeft = [subs[i].wallPointLeft for i = 1:N_subs]
        bottomWallNormal = [subs[i].bottomWallNormal for i = 1:N_subs]
        rightWallNormal = [subs[i].rightWallNormal for i = 1:N_subs]
        topWallNormal = [subs[i].topWallNormal for i = 1:N_subs]
        leftWallNormal = [subs[i].leftWallNormal for i = 1:N_subs]
        # return the new struct
        return new(N_subs,midPoints,Nx_coarse,Ny_coarse,
                    point1_coarse, point2_coarse, point3_coarse, point4_coarse, N_surfs_coarse, N_vols_coarse,
                    solidWalls, wallPointBottom, wallPointRight, wallPointTop, wallPointLeft,
                    bottomWallNormal, rightWallNormal, topWallNormal, leftWallNormal, Ndim, Ndim,
                    point1, point2, point3, point4, N_surfs, N_vols,
                    Area,Volume)        
    end
end
