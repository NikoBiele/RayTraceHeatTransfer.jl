function writeMatricesToCSV(FSS,FSG,FGS,FGG,Nx,Ny,N_rays,beta,omega,kappa,sigma_s)
    # This function saves the results of a ray tracing run,
    # which gives four exchange factor matrices.
    # These matrices are saved as CSV-files where the
    # name gives information on the settings that were
    # used to statistically measure them.

    # Surface to surface exchange factor matrix.
    # If the medium is transparent (absorbtion and scattering is zero) only this matrix has values.
    # If the medium is transparent this matrix gives the view factors of the enclosure.
    # Write FSS to csv-file.
    FSS_fileName = "Nx_$Nx.Ny_$Ny.FSS_$N_rays.Beta_$beta.Omega_$omega.Kappa_$kappa.Sigmas_$sigma_s.csv"
    CSV.write(FSS_fileName, Tables.table(FSS), header=false)

    # Surface to gas exchange factor matrix.
    # write FSG to csv-file.
    FSG_fileName = "Nx_$Nx.Ny_$Ny.FSG_$N_rays.Beta_$beta.Omega_$omega.Kappa_$kappa.Sigmas_$sigma_s.csv"
    CSV.write(FSG_fileName, Tables.table(FSG), header=false)

    # Gas to surface exchange factor matrix.
    # write FGS to csv-file.
    FGS_fileName = "Nx_$Nx.Ny_$Ny.FGS_$N_rays.Beta_$beta.Omega_$omega.Kappa_$kappa.Sigmas_$sigma_s.csv"
    CSV.write(FGS_fileName, Tables.table(FGS), header=false)

    # Gas to gas exchange factor matrices.
    # write FSG to csv-file.
    FGG_fileName = "Nx_$Nx.Ny_$Ny.FGG_$N_rays.Beta_$beta.Omega_$omega.Kappa_$kappa.Sigmas_$sigma_s.csv"
    CSV.write(FGG_fileName, Tables.table(FGG), header=false)
    
end