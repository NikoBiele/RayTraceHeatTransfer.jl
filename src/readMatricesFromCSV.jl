function readMatricesFromCSV(Nx,Ny,N_rays,beta,omega,kappa,sigma_s)
    # this function reads 4 exchange factor matrices from CSV-files
    # these files can be generated with the 'writeMatricesToCSV'-function

    # read FSS matrix
    FSS_fileName = "Nx_$Nx.Ny_$Ny.FSS_$N_rays.Beta_$beta.Omega_$omega.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FSS = Matrix(CSV.read(FSS_fileName, DataFrame, header=false))
    # read FSG matrix
    FSG_fileName = "Nx_$Nx.Ny_$Ny.FSG_$N_rays.Beta_$beta.Omega_$omega.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FSG = Matrix(CSV.read(FSG_fileName, DataFrame, header=false))
    # read FGS matrix
    FGS_fileName = "Nx_$Nx.Ny_$Ny.FGS_$N_rays.Beta_$beta.Omega_$omega.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FGS = Matrix(CSV.read(FGS_fileName, DataFrame, header=false))
    # read FGG matrix
    FGG_fileName = "Nx_$Nx.Ny_$Ny.FGG_$N_rays.Beta_$beta.Omega_$omega.Kappa_$kappa.Sigmas_$sigma_s.csv"
    FGG = Matrix(CSV.read(FGG_fileName, DataFrame, header=false))

    return FSS, FSG, FGS, FGG
end