module RayTraceHeatTransferMakieExt

using Makie
using RayTraceHeatTransfer
import RayTraceHeatTransfer: plotMesh2D, plotMesh3D, plotField3D

include(joinpath(@__DIR__, "plot2D", "plotMesh2D.jl"))
include(joinpath(@__DIR__, "plot3D", "plotMesh3D.jl"))
include(joinpath(@__DIR__, "plot3D", "plotField3D.jl"))

end