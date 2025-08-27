#####################################################
### THESE TESTS ARE OUTDATED AND WILL BE REPLACED ###
#####################################################

# using RayTraceHeatTransfer
# using Test
# using Interpolations

# ### To test this code I generate a mesh and ray trace it.
# ### In the end I test that the results are as they are expected to be (close to analytical).
# ### This way of testing utilize most of the functions in this package.
# ### Testing this package is a challenge since it relies on random numbers which introduces random noise.

# ### DEFINE COMMON GAS PROPERTIES

#     # emissivity of walls is 1 during tracing,
#     # but can be set differently during exchange factor calculations
#     # these are set to correspond to Crosbie & Schrenker
#     sigma_s = 0.0 # scattering coefficient
#     kappa = 1.0 # absorption coefficient
#     gas1 = GasProperties(sigma_s,kappa)

# ### DEFINE THE FIRST TEST GEOMETRY (Sides of square perpendicular to coordinate system)

#     println("Starting calculations for square aligned with coordinate system:")

#     # now we need to know the Coordinates of all the points in the enclosure
#     # define a vector of SubEnclosures
#     subs = SubEnclosure[]

#     # make one SubEnclosure at a time and push it into the vector
#     sub1 = SubEnclosure([0.0, 0.0],[1.0, 0.0],[1.0, 1.0],[0.0, 1.0],true,true,true,true)
#     push!(subs, sub1)

#     # generate mesh
#     Ndim = 11
#     mesh1 = RayTracingMesh(subs,Ndim);

# ### START THE RAY TRACING

#     displayWhileTracing = false
#     # number of rays to trace from each zone
#     N_rays_tot = 10_000_000 # we need several rays to decrease random noise
#     N_rays = trunc(Int, N_rays_tot/(mesh1.N_vols+mesh1.N_surfs))

#     # Here I make the calculation run in parallel on all available threads
#     nthreads = Threads.nthreads()

#     @time begin
#         println("Starting ray tracing of $N_rays_tot ray bundles in total ($N_rays per emitter):")
#         FSS, FSG, FGS, FGG = sampleDomain(mesh1,gas1,N_rays,nthreads,displayWhileTracing)
#         println("Ray tracing finished, the total time elapsed was:")
#     end;

# ### define common boundary conditions (we test emission from all four walls against analytical!)

#     # define which wall temperatures are fixed
#     # we only need to set those which should be fixed
#     # those which are set to negative will be interpreted as unknown
#     Tw_hot = 1000.0
#     qw_in = zeros(mesh1.N_subs,4)

#     # set the emissivities
#     epsw_in = ones(mesh1.N_subs,4)

#     # gas initial temperatures # set temperatures on a sub enclosure basis
#     Tg_in = zeros(mesh1.N_subs) .- 1 # unknown gas temperature
#     qg_in = zeros(mesh1.N_subs) #  radiative equilibrium

#     # set the emissivities
#     # convergence criteria (which ever happens first)
#     maxIter = 50
#     relTol = 1e-3

# ### TEST BOTTOM WALL EMISSIONS
#     Tw_in = zeros(mesh1.N_subs,4)
#     Tw_in[1,1] = Tw_hot # bottom wall
#     println("Calculating steady state temperature distribution (bottom wall emitter).")
#     Tw, Tg, qw, qg = steadyState(mesh1,FSS,FSG,FGS,FGG,
#                                         epsw_in,gas1,Tw_in,Tg_in,qw_in,qg_in);
#     Tg_matrix = plotTemperatureField(mesh1,Tg); #,Tw); # optional wall temperatures
#     Tg_bottom_1 = dropdims(Tg_matrix, dims=1)

# ### TEST RIGHT WALL EMISSIONS
#     Tw_in = zeros(mesh1.N_subs,4)
#     Tw_in[1,2] = Tw_hot # right wall
#     println("Calculating steady state temperature distribution (right wall emitter).")
#     Tw, Tg, qw, qg = steadyState(mesh1,FSS,FSG,FGS,FGG,
#                                         epsw_in,gas1,Tw_in,Tg_in,qw_in,qg_in);
#     Tg_matrix = plotTemperatureField(mesh1,Tg); #,Tw); # optional wall temperatures
#     Tg_right_1 = dropdims(Tg_matrix, dims=1)

# ### TEST TOP WALL EMISSIONS
#     Tw_in = zeros(mesh1.N_subs,4)
#     Tw_in[1,3] = Tw_hot # top wall
#     println("Calculating steady state temperature distribution (top wall emitter).")
#     Tw, Tg, qw, qg = steadyState(mesh1,FSS,FSG,FGS,FGG,
#                                         epsw_in,gas1,Tw_in,Tg_in,qw_in,qg_in);
#     Tg_matrix = plotTemperatureField(mesh1,Tg); #,Tw); # optional wall temperatures
#     Tg_top_1 = dropdims(Tg_matrix, dims=1)

# ### TEST LEFT WALL EMISSIONS
#     Tw_in = zeros(mesh1.N_subs,4)
#     Tw_in[1,4] = Tw_hot # left wall
#     println("Calculating steady state temperature distribution (left wall emitter).")
#     Tw, Tg, qw, qg = steadyState(mesh1,FSS,FSG,FGS,FGG,
#                                         epsw_in,gas1,Tw_in,Tg_in,qw_in,qg_in);
#     Tg_matrix = plotTemperatureField(mesh1,Tg); #,Tw); # optional wall temperatures
#     Tg_left_1 = dropdims(Tg_matrix, dims=1)

# ### DEFINE THE SECOND TEST GEOMETRY (Sides of square at a 45° angle with coordinate system)

#     println("Starting calculations for square at a 45° angle with coordinate system:")

#     # now we need to know the Coordinates of all the points in the enclosure
#     # define a vector of SubEnclosures
#     subs = SubEnclosure[]

#     # make one SubEnclosure at a time and push it into the vector
#     sub1 = SubEnclosure([0.0, 0.0],[sqrt(0.5), sqrt(0.5)],[0.0, 2*sqrt(0.5)],[-sqrt(0.5), sqrt(0.5)],true,true,true,true);
#     push!(subs, sub1)

#     # generate mesh
#     Ndim = 11
#     mesh1 = RayTracingMesh(subs,Ndim);

# ### START THE RAY TRACING

#     displayWhileTracing = false
#     # number of rays to trace from each zone
#     N_rays_tot = 10_000_000 # we need several rays to decrease random noise
#     N_rays = trunc(Int, N_rays_tot/(mesh1.N_vols+mesh1.N_surfs))

#     # Here I make the calculation run in parallel on all available threads
#     nthreads = Threads.nthreads()

#     @time begin
#         println("Starting ray tracing of $N_rays_tot ray bundles in total ($N_rays per emitter):")
#         FSS, FSG, FGS, FGG = sampleDomain(mesh1,gas1,N_rays,nthreads,displayWhileTracing)
#         println("Ray tracing finished, the total time elapsed was:")
#     end;

# ### define common boundary conditions (we test emission from all four walls against analytical!)

#     # define which wall temperatures are fixed
#     # we only need to set those which should be fixed
#     # those which are set to negative will be interpreted as unknown
#     Tw_hot = 1000.0
#     qw_in = zeros(mesh1.N_subs,4)

#     # set the emissivities
#     epsw_in = ones(mesh1.N_subs,4)

#     # gas initial temperatures # set temperatures on a sub enclosure basis
#     Tg_in = zeros(mesh1.N_subs) .- 1 # unknown gas temperature
#     qg_in = zeros(mesh1.N_subs) #  radiative equilibrium

#     # set the emissivities
#     # convergence criteria (which ever happens first)
#     maxIter = 50
#     relTol = 1e-3

# ### TEST BOTTOM WALL EMISSIONS
#     Tw_in = zeros(mesh1.N_subs,4)
#     Tw_in[1,1] = Tw_hot # bottom wall
#     println("Calculating steady state temperature distribution (bottom wall emitter).")
#     Tw, Tg, qw, qg = steadyState(mesh1,FSS,FSG,FGS,FGG,
#                                         epsw_in,gas1,Tw_in,Tg_in,qw_in,qg_in);
#     Tg_matrix = plotTemperatureField(mesh1,Tg); #,Tw); # optional wall temperatures
#     Tg_bottom_2 = dropdims(Tg_matrix, dims=1)

# ### TEST RIGHT WALL EMISSIONS
#     Tw_in = zeros(mesh1.N_subs,4)
#     Tw_in[1,2] = Tw_hot # right wall
#     println("Calculating steady state temperature distribution (right wall emitter).")
#     Tw, Tg, qw, qg = steadyState(mesh1,FSS,FSG,FGS,FGG,
#                                         epsw_in,gas1,Tw_in,Tg_in,qw_in,qg_in);
#     Tg_matrix = plotTemperatureField(mesh1,Tg); #,Tw); # optional wall temperatures
#     Tg_right_2 = dropdims(Tg_matrix, dims=1)

# ### TEST TOP WALL EMISSIONS
#     Tw_in = zeros(mesh1.N_subs,4)
#     Tw_in[1,3] = Tw_hot # top wall
#     println("Calculating steady state temperature distribution (top wall emitter).")
#     Tw, Tg, qw, qg = steadyState(mesh1,FSS,FSG,FGS,FGG,
#                                         epsw_in,gas1,Tw_in,Tg_in,qw_in,qg_in);
#     Tg_matrix = plotTemperatureField(mesh1,Tg); #,Tw); # optional wall temperatures
#     Tg_top_2 = dropdims(Tg_matrix, dims=1)

# ### TEST LEFT WALL EMISSIONS
#     Tw_in = zeros(mesh1.N_subs,4)
#     Tw_in[1,4] = Tw_hot # left wall
#     println("Calculating steady state temperature distribution (left wall emitter).")
#     Tw, Tg, qw, qg = steadyState(mesh1,FSS,FSG,FGS,FGG,
#                                         epsw_in,gas1,Tw_in,Tg_in,qw_in,qg_in);
#     Tg_matrix = plotTemperatureField(mesh1,Tg); #,Tw); # optional wall temperatures
#     Tg_left_2 = dropdims(Tg_matrix, dims=1)

# ### compare to analytical result

#     # check the center source function (1m x 1m, albedo = 1 which is equal to kappa = 1 for radiative equilibrium)
#     # relative optical depth
#     println("Comparing with analytical solution to validate correctness.")
#     relativeTauz = [0.0, 0.00611, 0.02037, 0.04251, 0.07216, 0.10884, 0.15194, 0.20076, 0.25449,
#                     0.31225, 0.37309, 0.43602, 0.50000, 0.56398, 0.62691, 0.68775, 0.74551,
#                     0.79924, 0.84806, 0.89116, 0.92784, 0.95749, 0.97963, 0.99390, 1.00000]
#     # value of the source function
#     sourceFuncCenter = [0.6293, 0.6198, 0.6017, 0.5767, 0.5460, 0.5108, 0.4724, 0.4323, 0.3919,
#                         0.3525, 0.3153, 0.2810, 0.2500, 0.2224, 0.1981, 0.1768, 0.1584, 0.1424,
#                         0.1287, 0.1171, 0.1073, 0.0992, 0.0930, 0.0885, 0.0863]

#     ### create intepolation object from Crosbie and Schrenker (1982) results
#     itp_CrosbieSchrenker = linear_interpolation(relativeTauz, sourceFuncCenter)

#     ### ray tracing source function along the center line
#     tauz_rays = 0.0+1/(Ndim*2):1/Ndim:1.0-1/(Ndim*2)

#     # intepolate Crosbie and Schrenker results to same position as ray tracing
#     CrosbieSchrenker_interpolated = itp_CrosbieSchrenker(tauz_rays) 

#     # calculate the dimensionless source function perpendicular to incident radiation
#     # for all four walls separately (square perpendicular to coordinate system)
#     sourceFun_rays_bottom_1 = (Tg_bottom_1[trunc(Int,1+(Ndim-1)/2),:]./Tw_hot).^4
#     sourceFun_rays_right_1 = (Tg_right_1[end:-1:1,trunc(Int,1+(Ndim-1)/2)]./Tw_hot).^4
#     sourceFun_rays_top_1 = (Tg_top_1[trunc(Int,1+(Ndim-1)/2),end:-1:1]./Tw_hot).^4
#     sourceFun_rays_left_1 = (Tg_left_1[:,trunc(Int,1+(Ndim-1)/2)]./Tw_hot).^4

#     # calculate the dimensionless source function perpendicular to incident radiation
#     # for all four walls separately (square at 45° angle with coordinate system)
#     sourceFun_rays_bottom_2 = (Tg_bottom_2[trunc(Int,1+(Ndim-1)/2),:]./Tw_hot).^4
#     sourceFun_rays_right_2 = (Tg_right_2[end:-1:1,trunc(Int,1+(Ndim-1)/2)]./Tw_hot).^4
#     sourceFun_rays_top_2 = (Tg_top_2[trunc(Int,1+(Ndim-1)/2),end:-1:1]./Tw_hot).^4
#     sourceFun_rays_left_2 = (Tg_left_2[:,trunc(Int,1+(Ndim-1)/2)]./Tw_hot).^4

# ### TEST VIEWFACTORS AGAINST THE EXAMPLES IN THE ORIGINAL ARTICLE

#     # "An analytic expression for radiation view factor between two arbitrarily oriented planar polygons"
#     # by Arvind Narayanaswamy
#     resultPaper = zeros(7)
#     F_AB = zeros(7)

#     # EXAMPLE 1
#     POLY_A = [0.0 0.0 0.0
#             1.0 0.0 0.0
#             1.0 1.0 0.0
#             0.0 1.0 0.0]
#     POLY_B = [0.0 0.0 1.0
#             1.0 0.0 1.0
#             1.0 1.0 1.0
#             0.0 1.0 1.0]
#     resultPaper[1] = 0.199825
#     F_AB[1], F_BA, area_A, area_B = viewFactor(POLY_A, POLY_B)

#     # EXAMPLE 2
#     POLY_A = [0.0 0.0 0.0
#             1.0 0.0 0.0
#             1.0 1.0 0.0
#             0.0 1.0 0.0]
#     POLY_B = [0.0 0.0 10.0
#             1.0 0.0 10.0
#             1.0 1.0 10.0
#             0.0 1.0 10.0]
#     resultPaper[2] = 3.16206e-3
#     F_AB[2], F_BA, area_A, area_B = viewFactor(POLY_A, POLY_B)

#     # EXAMPLE 3
#     POLY_A = [0.0 0.0 0.0
#             1.0 0.0 0.0
#             1.0 1.0 0.0
#             0.0 1.0 0.0]
#     POLY_B = [0.0 0.0 0.0
#             0.0 1.0 0.0
#             0.0 1.0 1.0
#             0.0 0.0 1.0]
#     resultPaper[3] = 0.200044
#     F_AB[3], F_BA, area_A, area_B = viewFactor(POLY_A, POLY_B)

#     # EXAMPLE 4
#     POLY_A = [0.0 0.0 0.0
#             0.0 1.0 0.0
#             1.0 1.0 0.0]
#     POLY_B = [1.0 0.0 1.0
#             1.0 1.0 1.0
#             0.0 1.0 1.0]
#     resultPaper[4] = 0.099912
#     F_AB[4], F_BA, area_A, area_B = viewFactor(POLY_A, POLY_B)

#     # EXAMPLE 5a
#     POLY_A = [0.0 0.5 0.0
#             1.0 0.0 0.0
#             1.0 1.0 0.0
#             0.0 1.0 0.0]
#     POLY_B = [2.0 0.5 0.0
#             3.0 0.0 0.5
#             3.0 2.0 0.5
#             2.0 1.5 0.0]
#     resultPaper[5] = 4.44228e-3
#     F_AB[5], F_BA, area_A, area_B = viewFactor(POLY_A, POLY_B)

#     # EXAMPLE 5b
#     POLY_A = [0.0 0.0 0.0
#             0.5 0.0 0.0
#             1.0 1.0 0.0
#             0.0 1.0 0.0]
#     POLY_B = [2.0 0.5 0.0
#             3.0 0.0 0.5
#             3.0 2.0 0.5
#             2.0 1.5 0.0]
#     resultPaper[6] = 3.63699e-3
#     F_AB[6], F_BA, area_A, area_B = viewFactor(POLY_A, POLY_B)

#     # EXAMPLE 6
#     POLY_A = [0.0 0.0 0.0
#             1.0 0.0 0.0
#             1.0 1.0 0.0]
#     POLY_B = [2.0 2.0 2.0
#             4.0 4.0 4.0
#             2.0 3.0 3.0]
#     resultPaper[7] = 1.06866e-3
#     F_AB[7], F_BA, area_A, area_B = viewFactor(POLY_A, POLY_B)

# ### TEST VALIDITY OF RESULTS

# @testset "RayTraceHeatTransfer.jl" begin
#     # if we are too far off the analytical results something probably went wrong
#     # perpendicular square
#     @test isapprox(CrosbieSchrenker_interpolated,sourceFun_rays_bottom_1) atol=0.02
#     @test isapprox(CrosbieSchrenker_interpolated,sourceFun_rays_right_1) atol=0.02
#     @test isapprox(CrosbieSchrenker_interpolated,sourceFun_rays_top_1) atol=0.02
#     @test isapprox(CrosbieSchrenker_interpolated,sourceFun_rays_left_1) atol=0.02
#     # square at 45° angle
#     @test isapprox(CrosbieSchrenker_interpolated,sourceFun_rays_bottom_2) atol=0.02
#     @test isapprox(CrosbieSchrenker_interpolated,sourceFun_rays_right_2) atol=0.02
#     @test isapprox(CrosbieSchrenker_interpolated,sourceFun_rays_top_2) atol=0.02
#     @test isapprox(CrosbieSchrenker_interpolated,sourceFun_rays_left_2) atol=0.02
#     # test view factor results
#     @test isapprox(resultPaper, F_AB) atol=1e-5
# end