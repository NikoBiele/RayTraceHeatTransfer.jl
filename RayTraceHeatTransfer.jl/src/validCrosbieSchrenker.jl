function validCrosbieSchrenker(N_rays_tot,Ndim)
    # this function performs a comparison of the ray traced result
    # with the analytical 2D solution tabulated by Crosbie & Schrenker (1982)

    # ensure that the dimension is an odd number (for comparison to analytical)
    if iseven(Ndim)
        error("Ndim should be an odd number")
    end

    ### DEFINE THE GEOMETRY

    # now we need to know the Coordinates of all the points in the enclosure
    # set the height of the outer points (the y-coordinates)
    displayGeometry = false;
    yLayersHeight = [0.0, 1.0];
    N_subs = length(yLayersHeight)-1; # number of sub-enclosures
    xLayersWidth = zeros(2, length(yLayersHeight));
    xLayersWidth[:,1] = [0.0, 1.0];
    xLayersWidth[:,2] = [0.0, 1.0];

    # coarse geometry
    # define the number of coarse splits in each enclosure
    Nx_coarse = 3; # must be minimum 3
    Ny_coarse = 2; # must be minimum 2
    point1_coarse, point2_coarse, point3_coarse, point4_coarse, N_surfs_coarse, N_vols_coarse =
                            geometry(yLayersHeight,xLayersWidth,Ny_coarse,Nx_coarse,displayGeometry);

    # define the number of fine splits in each enclosure
    Nx_fine = Ndim # must be minimum 3
    Ny_fine = Ndim
    point1_fine, point2_fine, point3_fine, point4_fine, N_surfs_fine, N_vols_fine =
                            geometry(yLayersHeight,xLayersWidth,Ny_fine,Nx_fine,displayGeometry);

    ### DEFINE GAS PROPERTIES

    # emissivity of walls is 1 during tracing, but can be set differently during exchange factor calculations
    # these are set to correspond to Crosbie & Schrenker
    sigma_s = 0.0
    kappa = 1.0
    beta = sigma_s+kappa
    omega = sigma_s/beta
    
    ### START THE RAY TRACING

    displayWhileTracing = false
    # number of rays to trace from each zone
    N_rays = trunc(Int, N_rays_tot/(Nx_fine*Ny_fine*N_subs+2*Nx_fine+2*Ny_fine))

    # Here I make the calculation run in parallel on all available threads
    if displayWhileTracing
        nthreads = 1 # 
    else
        nthreads = Threads.nthreads()
    end

    ### SAMPLE SURFACES
    println("ray tracing surface emissions")
    FSS, FSG = sampleSurfaces(point1_coarse, point2_coarse, point3_coarse, point4_coarse, Ny_coarse, Nx_coarse,
                        N_surfs_fine,N_vols_fine,point1_fine, point2_fine, point3_fine, point4_fine, Ny_fine, Nx_fine,
                        beta,omega,N_rays,displayWhileTracing,nthreads,N_subs);
    
    ### SAMPLE VOLUMES
    println("ray tracing volume emissions")
    FGS, FGG = sampleVolumes(point1_coarse, point2_coarse,point3_coarse, point4_coarse, Ny_coarse, Nx_coarse,
                        N_surfs_fine,N_vols_fine,point1_fine, point2_fine,point3_fine, point4_fine, Ny_fine, Nx_fine,
                        beta,omega,N_rays,displayWhileTracing,nthreads,N_subs);

    ### CALCULATE AREA AND VOLUME OF EACH ZONE

    width = 1.0 # width of domain
    Area, Volume = calculateAreaVolume(Nx_fine,Ny_fine,N_subs,width,point1_fine,point2_fine,point3_fine,point4_fine)

    ### CALCULATE THE STEADY STATE TEMPERATURE DISTRIBUTION

    # define which wall temperatures are fixed
    fixWalls = Vector{Bool}(undef, 4*Ndim)
    fixWalls .= true # all are fixed
    Tw_init = zeros(4*Ndim)
    Tw_init[1:Ny_fine] .= 0.0
    Tw_init[Ny_fine+1:2*Ny_fine] .= 0.0
    Tw_init[2*Ny_fine+1:2*Ny_fine+Nx_fine] .= 1.0
    Tw_init[2*Ny_fine+Nx_fine+1:2*Ny_fine+2*Nx_fine] .= 0.0
    # gas initial temperatures
    Tg_init = zeros(Ndim^2) .+ 100.0 # not fixed
    # set the emissivities
    epsw_vec = ones(4*Ndim)
    # convergence criteria (which ever happens first)
    maxIter = 200
    relTol = 1e-3

    println("Calculating steady state temperature distribution")
    Tw, Tg, break_count, Grelabs = steadyStateRigorous(Nx_fine,Ny_fine,N_subs,Area,Volume,FSS,FSG,FGS,FGG,
                                                        fixWalls,epsw_vec,kappa,maxIter,relTol,
                                                        Tw_init,Tg_init)

    display(plot(Grelabs[3:end],title="Convergence history"))
    sleep(5.0)

    println("plotting temperature distribution in the gas")
    Tg_matrix = Array{Float64}(undef, Nx_fine, Ny_fine*N_subs)
    Tg_count = 0
    for i = 1:Nx_fine
        for j = 1:Ny_fine*N_subs
            Tg_count += 1
            Tg_matrix[i,j] = Tg[Tg_count]
        end
    end
    display(contourf(Tg_matrix',aspect_ratio=1.0))
    sleep(5.0)


    ### COMPARISON WITH ANALYTICAL SOLUTION (Crosbie & Schrenker, 1982)
    
    println("comparison with analytical solution")
    # check the center source function (1m x 1m, albedo = 1 which is equal to kappa = 1 for radiative equilibrium)
    # relative optical depth
    relativeTauz = [0.0, 0.00611, 0.02037, 0.04251, 0.07216, 0.10884, 0.15194, 0.20076, 0.25449,
                    0.31225, 0.37309, 0.43602, 0.50000, 0.56398, 0.62691, 0.68775, 0.74551,
                    0.79924, 0.84806, 0.89116, 0.92784, 0.95749, 0.97963, 0.99390, 1.00000]
    # value of the source function
    sourceFuncCenter = [0.6293, 0.6198, 0.6017, 0.5767, 0.5460, 0.5108, 0.4724, 0.4323, 0.3919,
                        0.3525, 0.3153, 0.2810, 0.2500, 0.2224, 0.1981, 0.1768, 0.1584, 0.1424,
                        0.1287, 0.1171, 0.1073, 0.0992, 0.0930, 0.0885, 0.0863]

    display(plot(relativeTauz, sourceFuncCenter,label="Analytical Solution, Crosbie & Schenker (1982)"))
    N_rays_tot_mio = trunc(Int, N_rays_tot/1e6)
    label = "Monte Carlo Ray Tracing, $N_rays_tot_mio million rays"
    display(scatter!(0.0+1/(Ndim*2):1/Ndim:1.0-1/(Ndim*2),Tg_matrix[trunc(Int,(Ndim-1)/2),:].^4,marker=:star,label=label))
    display(xlabel!("Position / m"))
    display(ylabel!("Dimensionless Source Function"))

end