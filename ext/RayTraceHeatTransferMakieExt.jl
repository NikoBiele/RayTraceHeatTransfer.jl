module RayTraceHeatTransferMakieExt

using Makie
using RayTraceHeatTransfer
import RayTraceHeatTransfer: plotMesh, plotField

include("plotMesh.jl")
include("plotField.jl")

end