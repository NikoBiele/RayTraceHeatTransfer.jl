### Installation

To install the package, use the following command inside the Julia REPL:

```julia
using Pkg
Pkg.add("RayTraceHeatTransfer")
```

To load the package, use the command

```julia
using RayTraceHeatTransfer
```

### Enabling multithreading on the CPU

When performing Monte Carlo ray tracing it is advantageous to use multithreading.

If using VSCode, enable CPU multithreading by setting the following in settings.json:

```julia
  "julia.NumThreads": 16,
```

Instead of 16, choose the number you prefer.

If you're using Jupyter notebook run:

```julia
using IJulia
IJulia.installkernel("Julia 16 Threads", env=Dict(
    "JULIA_NUM_THREADS" => "16",
))
```

To confirm that the number of available threads correspond to the desired number run:

```
julia> Threads.nthreads()
```