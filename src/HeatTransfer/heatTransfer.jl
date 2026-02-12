# heat transfer
# black body
include(joinpath(@__DIR__, "blackBody", "emitFracBlackBodySpectrum.jl"))
include(joinpath(@__DIR__, "blackBody", "emitFracBlackBodySpectrumDerivative.jl"))
include(joinpath(@__DIR__, "blackBody", "getBinsEmissionFractions.jl"))
include(joinpath(@__DIR__, "blackBody", "solveTemperatureNewtonRaphson.jl"))
# equilibrium
include(joinpath(@__DIR__, "equilibrium", "WorkspaceStructs.jl"))
include(joinpath(@__DIR__, "equilibrium", "chooseSpectralMatrixType.jl"))
include(joinpath(@__DIR__, "equilibrium", "buildSystemMatrices.jl"))
include(joinpath(@__DIR__, "equilibrium", "equilibriumGrey2D.jl"))
include(joinpath(@__DIR__, "equilibrium", "equilibriumSpectral2D.jl"))
include(joinpath(@__DIR__, "equilibrium", "equilibriumSurfacesGrey3D.jl"))
include(joinpath(@__DIR__, "equilibrium", "equilibriumSurfacesSpectral3D.jl"))
include(joinpath(@__DIR__, "equilibrium", "setupBoundaryConditions.jl"))
include(joinpath(@__DIR__, "equilibrium", "solveEquilibrium.jl"))
include(joinpath(@__DIR__, "equilibrium", "updateSpectralEmission.jl"))
include(joinpath(@__DIR__, "equilibrium", "updateTemperaturesSpectral.jl"))
# exchange factor smoothing
include(joinpath(@__DIR__, "exchangeFactorSmoothing", "smoothExchangeFactors.jl"))
# heat transfer 2D
include(joinpath(@__DIR__, "writeResults", "writeResultsToDomain3D.jl"))
include(joinpath(@__DIR__, "writeResults", "writeTemperaturesHeatSources.jl"))