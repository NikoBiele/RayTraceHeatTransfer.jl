# direct ray tracing 2D
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing", "Emitter.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing", "directRayTracing.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing", "isotropicScatter.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing", "prepareEmitters.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing", "sampleReflectionDirection.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing", "traceSingleRay.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing", "updateHeatSource.jl"))

# exchange factors
include(joinpath(@__DIR__, "RayTracing2D", "ExchangeFactors", "exchangeRayTracing.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "ExchangeFactors", "parallelRayTracing.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "ExchangeFactors", "smoothExchangeFactors.jl"))

# shared functions
include(joinpath(@__DIR__, "RayTracing2D", "Shared", "traceRay.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared", "createIndexMapping.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared", "distToSurface.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared", "emitSurfaceRay.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared", "emitVolumeRay.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared", "findFace.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared", "findNearestFace.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared", "getGlobalIndex.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared", "getMappings.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared", "lambertSample3D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared", "multiDispatchRayTrace.jl"))

# ray tracing 3D
# empty for now

# view factors 3D
include(joinpath(@__DIR__, "ViewFactor3D", "Cl.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "edgePairParameters.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "f.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "fparallel.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "imagLi_2.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "viewFactor.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "exchangeFactors3D.jl"))