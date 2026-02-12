# direct ray tracing 2D
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing2D", "Emitter2D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing2D", "directRayTracing.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing2D", "isotropicScatter2D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing2D", "prepareEmitters.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing2D", "sampleReflectionDirection2D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing2D", "traceSingleRay.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "DirectTracing2D", "updateHeatSource.jl"))

# exchange factors
include(joinpath(@__DIR__, "RayTracing2D", "ExchangeFactors2D", "exchangeRayTracing.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "ExchangeFactors2D", "parallelRayTracing.jl"))

# shared functions
include(joinpath(@__DIR__, "RayTracing2D", "Shared2D", "traceRay.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared2D", "createIndexMapping2D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared2D", "distToSurface2D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared2D", "emitSurfaceRay2D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared2D", "emitVolumeRay2D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared2D", "findFace2D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared2D", "findNearestFace2D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared2D", "getGlobalIndex2D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared2D", "lambertSample2D.jl"))
include(joinpath(@__DIR__, "RayTracing2D", "Shared2D", "multiDispatchRayTrace2D.jl"))

# ray tracing 3D
# empty for now

# view factors 3D
include(joinpath(@__DIR__, "ViewFactor3D", "Cl3D.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "edgePairParameters3D.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "f3D.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "fparallel3D.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "imagLi2_3D.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "viewFactor3D.jl"))
include(joinpath(@__DIR__, "ViewFactor3D", "enclosureViewFactors3D.jl"))