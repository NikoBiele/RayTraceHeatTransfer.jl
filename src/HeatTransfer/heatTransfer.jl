# 2D
include(joinpath(@__DIR__, "heatTransfer2D", "emitFracBlackBodyElementSpectrum.jl"))
include(joinpath(@__DIR__, "heatTransfer2D", "emitFracBlackBodyElementSpectrumDeriv.jl"))
include(joinpath(@__DIR__, "heatTransfer2D", "getBinsEmissionFractions.jl"))
include(joinpath(@__DIR__, "heatTransfer2D", "setupBoundaryConditions.jl"))
include(joinpath(@__DIR__, "heatTransfer2D", "setupWavelengthBands.jl"))
include(joinpath(@__DIR__, "heatTransfer2D", "solveTemperatureNewtonRaphson.jl"))
include(joinpath(@__DIR__, "heatTransfer2D", "steadyStateGrey2D.jl"))
include(joinpath(@__DIR__, "heatTransfer2D", "steadyStateSpectral2D.jl"))
include(joinpath(@__DIR__, "heatTransfer2D", "updateSpectralEmission.jl"))
include(joinpath(@__DIR__, "heatTransfer2D", "updateTemperaturesSpectral.jl"))
include(joinpath(@__DIR__, "heatTransfer2D", "writeResultsToMesh.jl"))
# 3D
include(joinpath(@__DIR__, "heatTransfer3D", "steadyState3D.jl"))