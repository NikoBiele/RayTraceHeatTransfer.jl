using RayTraceHeatTransfer
using Test
using Interpolations

### To test this code I generate a mesh a ray trace it
### In the end I test that the results are as they are expected to be
### This way of testing utilize most of the functions in this package
### Testing this package is a challenge since it relies on random numbers which introduces random noise

# generate the a simple mesh to validate against the analytical result
yLayersHeight = [0.0, 1.0];
xLayersWidth = zeros(2, length(yLayersHeight));
xLayersWidth[:,1] = [0.0, 1.0];
xLayersWidth[:,2] = [0.0, 1.0];

# define the number of fine splits in each enclosure
Ndim = 11
Nx = Ndim # Ndim # must be minimum 3
Ny = Ndim #Ndim # must be minimum 2

# generate mesh
mesh1 = TracingMesh(Nx,Ny,xLayersWidth,yLayersHeight);

# define gas properties, like Crosbie & Schrenker
sigma_s = 0.0 # scattering coefficient
kappa = 1.0 # absorption coefficient
gas1 = GasProperties(sigma_s,kappa)

displayWhileTracing = false
# number of rays to trace from each zone
N_rays_tot = 1_000_000
N_rays = trunc(Int, N_rays_tot/(mesh1.N_vols+mesh1.N_surfs))

# Here I make the calculation run in parallel on all available threads
if displayWhileTracing
    nthreads = 1 # 
else
    nthreads = Threads.nthreads()
end

@time begin
    # generate the exchange factor matrices
    println("Starting ray tracing of $N_rays_tot ray bundles in total ($N_rays per emitter):")
    FSS, FSG, FGS, FGG = sampleDomain(mesh1,gas1,N_rays,nthreads,displayWhileTracing);
    println("Ray tracing finished, the total time elapsed was:")
end

### define boundary conditions

Tw_hot = 1000.0 # boundary condition bottom wall
# define which wall temperatures are fixed
fixWalls = Vector{Bool}(undef, 2*mesh1.Ny*mesh1.N_subs+2*mesh1.Nx)
fixWalls .= true # all are fixed
Tw_init = zeros(2*mesh1.Ny*mesh1.N_subs+2*mesh1.Nx)
Tw_init[1:mesh1.Ny*mesh1.N_subs] .= 0.0 # left wall
Tw_init[mesh1.Ny*mesh1.N_subs+1:2*mesh1.Ny*mesh1.N_subs] .= 0.0 # right wall
Tw_init[2*mesh1.Ny*mesh1.N_subs+1:2*mesh1.Ny*mesh1.N_subs+mesh1.Nx] .= Tw_hot # bottom wall
Tw_init[2*mesh1.Ny*mesh1.N_subs+mesh1.Nx+1:2*mesh1.Ny*mesh1.N_subs+2*Nx] .= 0.0 # top wall
# gas initial temperatures
Tg_init = zeros(mesh1.Ny*mesh1.Nx*mesh1.N_subs) .+ 0.0 # not fixed
# set the emissivities
epsw_vec = ones(2*mesh1.Ny*mesh1.N_subs+2*mesh1.Nx)
# convergence criteria (which ever happens first)
maxIter = 100
relTol = 1e-3

# calculate steady state temperatures
Tw, Tg, Gw, Gg, iter_count, Grelabs = steadyStateRigorous(mesh1,FSS,FSG,FGS,FGG,
                                                fixWalls,epsw_vec,kappa,maxIter,relTol,
                                                Tw_init,Tg_init);


Tg_matrix = plotTemperatureField(mesh1,Tg); #,Tw); # optional wall temperatures

### compare to analytical result

# check the center source function (1m x 1m, albedo = 1 which is equal to kappa = 1 for radiative equilibrium)
# relative optical depth
println("Comparing with analytical solution to validate correctness.")
relativeTauz = [0.0, 0.00611, 0.02037, 0.04251, 0.07216, 0.10884, 0.15194, 0.20076, 0.25449,
                0.31225, 0.37309, 0.43602, 0.50000, 0.56398, 0.62691, 0.68775, 0.74551,
                0.79924, 0.84806, 0.89116, 0.92784, 0.95749, 0.97963, 0.99390, 1.00000]
# value of the source function
sourceFuncCenter = [0.6293, 0.6198, 0.6017, 0.5767, 0.5460, 0.5108, 0.4724, 0.4323, 0.3919,
                    0.3525, 0.3153, 0.2810, 0.2500, 0.2224, 0.1981, 0.1768, 0.1584, 0.1424,
                    0.1287, 0.1171, 0.1073, 0.0992, 0.0930, 0.0885, 0.0863]

### create intepolation object from Crosbie and Schrenker (1982) results
itp_CrosbieSchrenker = linear_interpolation(relativeTauz, sourceFuncCenter)

### ray tracing source function along the center line
tauz_rays = 0.0+1/(Ndim*2):1/Ndim:1.0-1/(Ndim*2)
sourceFun_rays = (Tg_matrix[trunc(Int,(Ndim-1)/2),:]./Tw_hot).^4

# intepolate Crosbie and Schrenker results to same position as ray tracing
CrosbieSchrenker_interpolated = itp_CrosbieSchrenker(tauz_rays) 

@testset "RayTraceHeatTransfer.jl" begin
    # if we are too far off the analytical results something probably went wrong
    @test isapprox(CrosbieSchrenker_interpolated,sourceFun_rays) atol=0.05
end
