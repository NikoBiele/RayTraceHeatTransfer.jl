# direct ray tracing 2D
include("RayTracing2D/DirectTracing/Emitter.jl")
include("RayTracing2D/DirectTracing/directRayTracing.jl")
include("RayTracing2D/DirectTracing/isotropicScatter.jl")
include("RayTracing2D/DirectTracing/prepareEmitters.jl")
include("RayTracing2D/DirectTracing/sampleReflectionDirection.jl")
include("RayTracing2D/DirectTracing/traceSingleRay.jl")
include("RayTracing2D/DirectTracing/updateHeatSource.jl")

# exchange factors
include("RayTracing2D/ExchangeFactors/checkEnergyConservation.jl")
include("RayTracing2D/ExchangeFactors/exchangeRayTracing.jl")
include("RayTracing2D/ExchangeFactors/parallelRayTracing.jl")
include("RayTracing2D/ExchangeFactors/smoothExchangeFactors.jl")
include("RayTracing2D/ExchangeFactors/checkReciprocity.jl")

# shared functions
include("RayTracing2D/Shared/traceRay.jl")
include("RayTracing2D/Shared/createIndexMapping.jl")
include("RayTracing2D/Shared/distToSurface.jl")
include("RayTracing2D/Shared/emitSurfaceRay.jl")
include("RayTracing2D/Shared/emitVolumeRay.jl")
include("RayTracing2D/Shared/findFace.jl")
include("RayTracing2D/Shared/findNearestFace.jl")
include("RayTracing2D/Shared/getGlobalIndex.jl")
include("RayTracing2D/Shared/getMappings.jl")
include("RayTracing2D/Shared/lambertSample3D.jl")
include("RayTracing2D/Shared/multiDispatchRayTrace.jl")

# ray tracing 3D
# empty for now

# view factors 3D
include("ViewFactor3D/Cl.jl")
include("ViewFactor3D/edgePairParameters.jl")
include("ViewFactor3D/f.jl")
include("ViewFactor3D/fparallel.jl")
include("ViewFactor3D/imagLi_2.jl")
include("ViewFactor3D/viewFactor.jl")
include("ViewFactor3D/exchangeFactors3D.jl")