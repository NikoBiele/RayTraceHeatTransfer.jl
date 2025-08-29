module RayTraceHeatTransferPlotsExt

using Plots
using RayTraceHeatTransfer
import RayTraceHeatTransfer: plotField2D

include(joinpath(@__DIR__, "plot2D", "plotField2D.jl"))

end