struct TracingMesh

    # number of sub-enclosures
    N_subs::Int64

    # the coarse mesh
    Nx_coarse::Int64
    Ny_coarse::Int64
    point1_coarse::Matrix{SVector{2,Float64}}
    point2_coarse::Matrix{SVector{2,Float64}}
    point3_coarse::Matrix{SVector{2,Float64}}
    point4_coarse::Matrix{SVector{2,Float64}}
    N_surfs_coarse::Int64
    N_vols_coarse::Int64
    NeighborIndices_coarse::Matrix{Array{SVector{2, Int64}}}
    solidWalls::Vector{Vector{Bool}}

    # the fine mesh
    Nx::Int64
    Ny::Int64
    point1::Matrix{SVector{2,Float64}}
    point2::Matrix{SVector{2,Float64}}
    point3::Matrix{SVector{2,Float64}}
    point4::Matrix{SVector{2,Float64}}
    N_surfs::Int64
    N_vols::Int64
    NeighborIndices::Matrix{Array{SVector{2, Int64}}}
    Area::Vector{Float64}
    Volume::Vector{Float64}

    # generate the fields automatically from fewer inputs about the bounding geometry
    function TracingMesh(Nx::Int64,Ny::Int64,xLayersWidth::Matrix{Float64},yLayersHeight::Vector{Float64})
        # number of sub-enclosures
        N_subs = length(yLayersHeight)-1
        # define the number of coarse splits in each enclosure
        # these are not allowed to be user defined, since they are very performance critical
        Nx_coarse = 1
        Ny_coarse = 1
        # mesh to get coarse geometry
        point1_coarse, point2_coarse, point3_coarse, point4_coarse, N_surfs_coarse, N_vols_coarse, NeighborIndices_coarse =
                            meshGeometry(yLayersHeight,xLayersWidth,Ny_coarse,Nx_coarse);
        # mesh to get fine geometry
        point1, point2, point3, point4, N_surfs, N_vols, NeighborIndices =
                            meshGeometry(yLayersHeight,xLayersWidth,Ny,Nx);
        # define the outer solid walls
        solidWalls = defineSolidWalls(yLayersHeight::Vector{Float64})
        # calculate the area and volume of each cell in the fine mesh
        Area, Volume = calculateAreaVolume(Nx,Ny,N_subs,point1,point2,point3,point4);
        # return the new struct
        return new(N_subs,Nx_coarse,Ny_coarse,
                    point1_coarse, point2_coarse, point3_coarse, point4_coarse, N_surfs_coarse, N_vols_coarse,
                    NeighborIndices_coarse, solidWalls,
                    Nx,Ny,
                    point1, point2, point3, point4, N_surfs, N_vols, NeighborIndices,
                    Area,Volume)        
    end

end

struct GasProperties
    sigma_s::Float64
    kappa::Float64
    beta::Float64
    omega::Float64
    function GasProperties(sigma_s::Float64,kappa::Float64)
        beta = sigma_s+kappa
        omega = sigma_s/beta 
        return new(sigma_s,kappa,beta,omega)
    end
end