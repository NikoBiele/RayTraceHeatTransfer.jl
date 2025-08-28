# 2D
include(joinpath(@__DIR__, "Domain2D", "domain2D.jl"))
include(joinpath(@__DIR__, "Domain2D", "meshing2D.jl"))
# 3D
include(joinpath(@__DIR__, "Domain3D_faces", "domain3D_faces.jl"))
include(joinpath(@__DIR__, "Domain3D_faces", "meshing3D.jl"))
include(joinpath(@__DIR__, "Domain3D_faces", "projectPlane.jl"))