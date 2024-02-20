function validCrosbieSchrenker(N_rays_tot,Ndim,Tw_hot)
    # this function performs a comparison of the ray traced result
    # with the analytical 2D solution tabulated by Crosbie & Schrenker (1982)

    # ensure that the dimension is an odd number (for comparison to analytical)
    if iseven(Ndim)
        error("Ndim should be an odd number to validate against the analytical result.")
    end

    ### DEFINE THE GEOMETRY

    # now we need to know the Coordinates of all the points in the enclosure
    # set the height of the outer points (the y-coordinates)
    yLayersHeight = [0.0, 1.0];
    xLayersWidth = zeros(2, length(yLayersHeight));
    xLayersWidth[:,1] = [0.0, 1.0];
    xLayersWidth[:,2] = [0.0, 1.0];

    # define the number of fine splits in each enclosure
    Nx = Ndim
    Ny = Ndim

    # generate geometry
    mesh1 = TracingMesh(Nx,Ny,xLayersWidth,yLayersHeight);

    # plot the geometry
    displayMesh(mesh1)
    println("Waiting 5 seconds to view mesh...")
    sleep(5.0)
    
    ### DEFINE GAS PROPERTIES

    # emissivity of walls is 1 during tracing, but can be set differently during exchange factor calculations
    # these are set to correspond to Crosbie & Schrenker
    sigma_s = 0.0 # scattering coefficient
    kappa = 1.0 # absorption coefficient
    gas1 = GasProperties(sigma_s,kappa)
    
    ### START THE RAY TRACING

    displayWhileTracing = false
    # number of rays to trace from each zone
    N_rays = trunc(Int, N_rays_tot/(mesh1.N_vols+mesh1.N_surfs))

    # Here I make the calculation run in parallel on all available threads
    if displayWhileTracing
        nthreads = 1 # 
    else
        nthreads = Threads.nthreads()
    end

    @time begin
        println("Starting ray tracing of $N_rays_tot ray bundles in total ($N_rays per emitter):")
        FSS, FSG, FGS, FGG = sampleDomain(mesh1,gas1,N_rays,nthreads,displayWhileTracing)
        println("Ray tracing finished, the total time elapsed was:")
    end;

    ### CALCULATE THE STEADY STATE TEMPERATURE DISTRIBUTION

    # define which wall temperatures are fixed
    fixWalls = Vector{Bool}(undef, 2*mesh1.Ny*mesh1.N_subs+2*mesh1.Nx)
    fixWalls .= true # all are fixed
    Tw_init = zeros(2*mesh1.Ny*mesh1.N_subs+2*mesh1.Nx)
    Tw_init[1:mesh1.Ny*mesh1.N_subs] .= 0.0 # left wall
    Tw_init[mesh1.Ny*mesh1.N_subs+1:2*mesh1.Ny*mesh1.N_subs] .= 0.0 # right wall
    Tw_init[2*mesh1.Ny*mesh1.N_subs+1:2*mesh1.Ny*mesh1.N_subs+mesh1.Nx] .= Tw_hot # bottom wall
    Tw_init[2*mesh1.Ny*mesh1.N_subs+mesh1.Nx+1:2*mesh1.Ny*mesh1.N_subs+2*Nx] .= 0.0 # top wall
    # gas initial temperatures
    Tg_init = zeros(mesh1.Ny*mesh1.Nx*mesh1.N_subs) .+ 100.0 # not fixed
    # set the emissivities
    epsw_vec = ones(2*mesh1.Ny*mesh1.N_subs+2*mesh1.Nx)
    # convergence criteria (which ever happens first)
    maxIter = 100
    relTol = 1e-3

    println("Calculating steady state temperature distribution.")
    Tw, Tg, Gw, Gg, iter_count, Grelabs = steadyStateRigorous(mesh1,FSS,FSG,FGS,FGG,
                                                        fixWalls,epsw_vec,kappa,maxIter,relTol,
                                                        Tw_init,Tg_init);

    display(plot(Grelabs[3:end],title="Convergence history",legend=false))
    display(xlabel!("Iterations"))
    display(ylabel!("Relative change in heat flux"))
    println("Waiting 5 seconds to view convergence history...")
    sleep(5.0)

    println("Plotting temperature distribution in the gas.")
    Tg_matrix = plotTemperatureField(mesh1,Tg); #,Tw); # optional wall temperatures
    println("Waiting 5 seconds to view temperature distribution...")
    sleep(5.0)

    ### COMPARISON WITH ANALYTICAL SOLUTION (Crosbie & Schrenker, 1982)
    
    println("Comparison with analytical solution.")
    # check the center source function (1m x 1m, albedo = 1 which is equal to kappa = 1 for radiative equilibrium)
    # relative optical depth
    relativeTauz = [0.0, 0.00611, 0.02037, 0.04251, 0.07216, 0.10884, 0.15194, 0.20076, 0.25449,
                    0.31225, 0.37309, 0.43602, 0.50000, 0.56398, 0.62691, 0.68775, 0.74551,
                    0.79924, 0.84806, 0.89116, 0.92784, 0.95749, 0.97963, 0.99390, 1.00000]
    # value of the source function
    sourceFuncCenter = [0.6293, 0.6198, 0.6017, 0.5767, 0.5460, 0.5108, 0.4724, 0.4323, 0.3919,
                        0.3525, 0.3153, 0.2810, 0.2500, 0.2224, 0.1981, 0.1768, 0.1584, 0.1424,
                        0.1287, 0.1171, 0.1073, 0.0992, 0.0930, 0.0885, 0.0863]

    display(plot(relativeTauz, sourceFuncCenter,label="Analytical Solution, Crosbie & Schrenker (1982)"))
    N_rays_tot_mio = trunc(Int, N_rays_tot/1e6)
    label = "Monte Carlo Ray Tracing, $N_rays_tot_mio million rays"
    display(scatter!(0.0+1/(Ndim*2):1/Ndim:1.0-1/(Ndim*2),(Tg_matrix[trunc(Int,(Ndim-1)/2),:]./Tw_hot).^4,marker=:star,label=label))
    display(xlabel!("Position / m"))
    display(ylabel!("Dimensionless Source Function"))

end