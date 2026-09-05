module RayTraceHeatTransferMakieExt

using GeometryBasics
using Makie
using RayTraceHeatTransfer
import RayTraceHeatTransfer: SurfaceDomain3D, plotDomain3D, plotMesh, DomainPlot3D

include("plotMesh.jl")
include("plotDomain3D.jl")

end