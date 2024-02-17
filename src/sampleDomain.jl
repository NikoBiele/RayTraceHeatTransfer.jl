function sampleDomain(point1_coarse::Matrix{SVector{2,Float64}}, point2_coarse::Matrix{SVector{2,Float64}},
                        point3_coarse::Matrix{SVector{2,Float64}}, point4_coarse::Matrix{SVector{2,Float64}}, Ny_coarse::Int, Nx_coarse::Int,
                        N_surfs_fine::Int,N_vols_fine::Int,point1_fine::Matrix{SVector{2,Float64}}, point2_fine::Matrix{SVector{2,Float64}},
                        point3_fine::Matrix{SVector{2,Float64}}, point4_fine::Matrix{SVector{2,Float64}}, Ny_fine::Int, Nx_fine::Int,
                        beta::Float64,omega::Float64,N_rays::Int64,displayWhileTracing::Bool,nthreads::Int,N_subs::Int,
                        NeighborIndices_coarse::Matrix{Array{SVector{2, Int64}}})

    ### SAMPLE SURFACES
    FSS, FSG = sampleSurfaces(point1_coarse,point2_coarse,point3_coarse,point4_coarse,Ny_coarse,Nx_coarse,
                            N_surfs_fine,N_vols_fine,point1_fine,point2_fine,point3_fine,point4_fine,Ny_fine,Nx_fine,
                            beta,omega,N_rays,displayWhileTracing,nthreads,N_subs,NeighborIndices_coarse);

    ### SAMPLE VOLUMES
    FGS, FGG = sampleVolumes(point1_coarse,point2_coarse,point3_coarse,point4_coarse,Ny_coarse,Nx_coarse,
                            N_surfs_fine,N_vols_fine,point1_fine,point2_fine,point3_fine,point4_fine,Ny_fine,Nx_fine,
                            beta,omega,N_rays,displayWhileTracing,nthreads,N_subs,NeighborIndices_coarse);

    return FSS, FSG, FGS, FGG
end