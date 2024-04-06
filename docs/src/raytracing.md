### Ray trace the geometry

In the previous section the geometry was created.
Now it is time to sample the geometry.
This package works by sampling every element of the entire geometry a specified number of times.
This sampling process generates 'exchange factor' matrices.
The exchange factors describe the connectivity of the domain.
This connectivity depends on the properties of the participating medium (the enclosed gas).

The first thing we need to do is to specify the properties of the participating medium.
This is done by speciying the scattering coefficient and the absorption coefficient and generating an instance of the GasProperties type.

```julia
sigma_s = 0.0 # scattering coefficient
kappa = 1.0 # absorption coefficient
gas1 = GasProperties(sigma_s,kappa)
```

Then we specify the total number of rays that should be traced in the domain.
More rays give a more accurate result but is also more computationally demanding.
Here we choose 10 million rays.

```julia
N_rays_tot = 10_000_000
```

Then we divide these rays among the total number of elements in the domain

```julia
N_rays = trunc(Int, N_rays_tot/(mesh1.N_vols+mesh1.N_surfs))
```

The sampling process gives the user the option of viewing the ray tracing 'live' as it takes place.
This option is however very computationally demanding and must execute single threaded.
Therefore, if enabling this option, only a few rays should be sampled from each emitter.
But we will leave this option turned off.

```julia
displayWhileTracing = false
```

Next we make the calculation run in parallel on all available threads.
This is done to accelerate the ray tracing.
The ray tracing can run in parallel on any number of threads.

```julia
if displayWhileTracing
    nthreads = 1
else
    nthreads = Threads.nthreads()
end
```

Lastly we sample the domain and time the entire process.
This sampling process generates the exchange factor matrices.
The sampling function prints its progress to the Julia REPL.

```julia
@time begin
    println("Starting ray tracing of $N_rays_tot ray bundles in total ($N_rays per emitter):")
    FSS, FSG, FGS, FGG = sampleDomain(mesh1,gas1,N_rays,nthreads,displayWhileTracing)
    println("Ray tracing finished, the total time elapsed was:")
end;
```

It is recommended (but not mandatory) to save the exchange factor matrices to disk before continuing, to prevent loss of the data from a computer crash.
This can be achieved in the following way (matrices are written to csv-files in the current working directory):

```julia
writeMatricesToCSV(FSS,FSG,FGS,FGG,mesh1,gas1,N_rays)
```

Now that the domain has been sampled and the exchange factors has been statistically measured the next step is to solve heat transfer problems.