# heat transfer
# black body
include(joinpath(@__DIR__, "blackBody", "emitFracBlackBodySpectrum.jl"))
include(joinpath(@__DIR__, "blackBody", "emitFracBlackBodySpectrumDerivative.jl"))
include(joinpath(@__DIR__, "blackBody", "getBinsEmissionFractions.jl"))
include(joinpath(@__DIR__, "blackBody", "solveTemperatureNewtonRaphson.jl"))
include(joinpath(@__DIR__, "blackBody", "planckTable.jl"))
# equilibrium
include(joinpath(@__DIR__, "equilibrium", "WorkspaceStructs.jl"))
include(joinpath(@__DIR__, "equilibrium", "buildSystemMatrix.jl"))
include(joinpath(@__DIR__, "equilibrium", "equilibriumGrey2D.jl"))
include(joinpath(@__DIR__, "equilibrium", "equilibriumSpectral2D.jl"))
include(joinpath(@__DIR__, "equilibrium", "equilibriumSurfacesGrey3D.jl"))
include(joinpath(@__DIR__, "equilibrium", "equilibriumSurfacesSpectral3D.jl"))
include(joinpath(@__DIR__, "equilibrium", "setupBoundaryConditions.jl"))
include(joinpath(@__DIR__, "equilibrium", "solveEquilibrium.jl"))
include(joinpath(@__DIR__, "equilibrium", "updateSpectralEmission.jl"))
include(joinpath(@__DIR__, "equilibrium", "updateTemperaturesSpectral.jl"))
# heat transfer 2D
include(joinpath(@__DIR__, "writeResults", "writeResultsToDomain3D.jl"))
include(joinpath(@__DIR__, "writeResults", "writeTemperaturesHeatSources.jl"))
# adaptive spectral binning
include(joinpath(@__DIR__, "spectral", "SpectralModels.jl"))
include(joinpath(@__DIR__, "spectral", "spectralHooks.jl"))
include(joinpath(@__DIR__, "spectral", "adaptiveSpectralBins.jl"))
# directional (angular) models
include(joinpath(@__DIR__, "directional", "DirectionalModels.jl"))
include(joinpath(@__DIR__, "equilibrium", "equilibriumGreyDirectional2D.jl"))
include(joinpath(@__DIR__, "equilibrium", "equilibriumSpectralDirectional2D.jl"))