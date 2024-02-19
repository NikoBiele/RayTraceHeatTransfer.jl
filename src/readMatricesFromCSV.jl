function readMatricesFromCSV(N_subs,Nx,Ny,kappa,sigma_s,N_rays)
    # this function reads 4 exchange factor matrices from CSV-files
    # these files can be generated with the 'writeMatricesToCSV'-function

    # read FSS matrix
    FSS_fileName = "FSS_$N_rays.N_subs_$N_subs.Nx_$Nx.Ny_$Ny.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FSS = Matrix(CSV.read(FSS_fileName, DataFrame, header=false))
    # read FSG matrix
    FSG_fileName = "FSG_$N_rays.N_subs_$N_subs.Nx_$Nx.Ny_$Ny.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FSG = Matrix(CSV.read(FSG_fileName, DataFrame, header=false))
    # read FGS matrix
    FGS_fileName = "FGS_$N_rays.N_subs_$N_subs.Nx_$Nx.Ny_$Ny.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FGS = Matrix(CSV.read(FGS_fileName, DataFrame, header=false))
    # read FGG matrix
    FGG_fileName = "FGG_$N_rays.N_subs_$N_subs.Nx_$Nx.Ny_$Ny.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FGG = Matrix(CSV.read(FGG_fileName, DataFrame, header=false))

    return FSS, FSG, FGS, FGG
end