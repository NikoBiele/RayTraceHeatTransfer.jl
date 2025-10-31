# 2D
include(joinpath(@__DIR__, "RayTracingDomain2D", "RayTracingDomain2D.jl"))
include(joinpath(@__DIR__, "RayTracingDomain2D", "meshing2D.jl"))
# 3D
include(joinpath(@__DIR__, "ViewFactorDomain3D", "ViewFactorDomain3D.jl"))
include(joinpath(@__DIR__, "ViewFactorDomain3D", "meshing3D.jl"))
include(joinpath(@__DIR__, "ViewFactorDomain3D", "projectPlane.jl"))