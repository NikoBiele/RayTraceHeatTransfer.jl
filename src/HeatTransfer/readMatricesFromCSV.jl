"""
    readMatricesFromCSV(N_subs::Int64,Ndim::Int64,kappa::Float64,
                            sigma_s::Float64,N_rays::Int64)

This function reads 4 exchange factor matrices from CSV-files.
These files can be generated with the 'writeMatricesToCSV'-function.
"""
function readMatricesFromCSV(N_subs::Int64,Ndim::Int64,kappa::Float64,
                            sigma_s::Float64,N_rays::Int64)

    # read FSS matrix
    FSS_fileName = "FSS_$N_rays.N_subs_$N_subs.Ndim_$Ndim.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FSS = Matrix(CSV.read(FSS_fileName, DataFrame, header=false))
    # read FSG matrix
    FSG_fileName = "FSG_$N_rays.N_subs_$N_subs.Ndim_$Ndim.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FSG = Matrix(CSV.read(FSG_fileName, DataFrame, header=false))
    # read FGS matrix
    FGS_fileName = "FGS_$N_rays.N_subs_$N_subs.Ndim_$Ndim.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FGS = Matrix(CSV.read(FGS_fileName, DataFrame, header=false))
    # read FGG matrix
    FGG_fileName = "FGG_$N_rays.N_subs_$N_subs.Ndim_$Ndim.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FGG = Matrix(CSV.read(FGG_fileName, DataFrame, header=false))

    return FSS, FSG, FGS, FGG
end