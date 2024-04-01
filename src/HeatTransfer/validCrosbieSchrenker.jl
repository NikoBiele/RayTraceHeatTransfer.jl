"""
    validCrosbieSchrenker(N_rays_tot::Int64,Ndim::Int64,Tw_hot::Float64)

This function performs a comparison of the ray tracing result in a square
with the analytical 2D solution tabulated by Crosbie & Schrenker (1982).
"""
function validCrosbieSchrenker(N_rays_tot::Int64,Ndim::Int64,Tw_hot::Float64)
    
    # ensure that the dimension is an odd number (for comparison to analytical)
    if iseven(Ndim)
        error("Ndim should be an odd number to validate against the analytical result.")
    end

    ### DEFINE THE GEOMETRY

    # now we need to know the Coordinates of all the points in the enclosure
    # define a vector of SubEnclosures
    subs = SubEnclosure[]

    # make one SubEnclosure at a time and push it into the vector
    sub1 = SubEnclosure([0.0, 0.0],[1.0, 0.0],[1.0, 1.0],[0.0, 1.0],true,true,true,true)
    push!(subs, sub1)

    # generate mesh
    mesh1 = RayTracingMesh(subs,Ndim);

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

    # we change so that temperatures and source terms
    # are set on a sub-enclosure basis
    # then we just have to 'unpack' these within the function

    # define which wall temperatures are fixed
    # we only need to set those which should be fixed
    # those which are set to negative will be interpreted as unknown
    Tw_in = zeros(mesh1.N_subs,4)
    qw_in = zeros(mesh1.N_subs,4)
    Tw_in[1,mesh1.solidWalls[1]] .= 0.0
    Tw_in[1,1] = Tw_hot
    # qw_in[1,mesh1.solidWalls[1]] .= 0.0 # unknown wall flux

    # set the emissivities
    epsw_in = ones(mesh1.N_subs,4)

    # gas initial temperatures # set temperatures on a sub enclosure basis
    Tg_in = zeros(mesh1.N_subs) .- 1 # unknown gas temperature
    qg_in = zeros(mesh1.N_subs) #  radiative equilibrium

    # set the emissivities
    # convergence criteria (which ever happens first)
    maxIter = 50
    relTol = 1e-3

    println("Calculating steady state temperature distribution.")
    Tw, Tg, Gw, Gg, iter_count, Grelabs = steadyStateRigorous(mesh1,FSS,FSG,FGS,FGG,
                                                        epsw_in,gas1,maxIter,relTol,
                                                        Tw_in,Tg_in,qw_in,qg_in);

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
    Tg_matrix = dropdims(Tg_matrix, dims=1)
    display(scatter!(0.0+1/(Ndim*2):1/Ndim:1.0-1/(Ndim*2),(Tg_matrix[trunc(Int,1+(Ndim-1)/2),:]./Tw_hot).^4,marker=:star,label=label))
    display(xlabel!("Position / m"))
    display(ylabel!("Dimensionless Source Function"))

end