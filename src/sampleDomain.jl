function sampleDomain(mesh::TracingMesh,gas::GasProperties,N_rays::Int64,nthreads::Int64,displayWhileTracing::Bool)
    # point1_coarse::Matrix{SVector{2,Float64}}, point2_coarse::Matrix{SVector{2,Float64}},
    #                     point3_coarse::Matrix{SVector{2,Float64}}, point4_coarse::Matrix{SVector{2,Float64}}, Ny_coarse::Int, Nx_coarse::Int,
    #                     N_surfs::Int,N_vols::Int,point1::Matrix{SVector{2,Float64}}, point2::Matrix{SVector{2,Float64}},
    #                     point3::Matrix{SVector{2,Float64}}, point4::Matrix{SVector{2,Float64}}, Ny::Int, Nx::Int,
    #                     beta::Float64,omega::Float64,N_rays::Int64,displayWhileTracing::Bool,nthreads::Int,N_subs::Int,
    #                     NeighborIndices_coarse::Matrix{Array{SVector{2, Int64}}})

    ### SAMPLE SURFACES
    FSS, FSG = sampleSurfaces(mesh,gas,N_rays,nthreads,displayWhileTracing);

    ### SAMPLE VOLUMES
    FGS, FGG = sampleVolumes(mesh,gas,N_rays,nthreads,displayWhileTracing);

    return FSS, FSG, FGS, FGG
end