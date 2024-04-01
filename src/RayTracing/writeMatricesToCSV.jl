"""
    writeMatricesToCSV(FSS::Matrix{Float64},FSG::Matrix{Float64},
                        FGS::Matrix{Float64},FGG::Matrix{Float64},
                        mesh::RayTracingMesh,gas::GasProperties,
                        N_rays::Int64)

This function saves the results of a ray tracing run,
which gives four exchange factor matrices.
These matrices are saved as CSV-files where the
name gives information on the settings that were
used to statistically measure them.
"""
function writeMatricesToCSV(FSS::Matrix{Float64},FSG::Matrix{Float64},
                            FGS::Matrix{Float64},FGG::Matrix{Float64},
                            mesh::RayTracingMesh,gas::GasProperties,
                            N_rays::Int64)
    # extract data for naming of the files
    N_subs = mesh.N_subs
    Ndim = mesh.Nx
    kappa = gas.kappa
    sigma_s = gas.sigma_s

    # Surface to surface exchange factor matrix.
    # If the medium is transparent (absorbtion and scattering is zero) only this matrix has values.
    # If the medium is transparent this matrix gives the view factors of the enclosure (almost, since we have diffuse emission).
    # Write FSS to csv-file.
    FSS_fileName = "FSS_$N_rays.N_subs_$N_subs.Ndim_$Ndim.Kappa_$kappa.Sigmas_$sigma_s.csv"
    CSV.write(FSS_fileName, Tables.table(FSS), header=false)

    # Surface to gas exchange factor matrix.
    # write FSG to csv-file.
    FSG_fileName = "FSG_$N_rays.N_subs_$N_subs.Ndim_$Ndim.Kappa_$kappa.Sigmas_$sigma_s.csv"
    CSV.write(FSG_fileName, Tables.table(FSG), header=false)

    # Gas to surface exchange factor matrix.
    # write FGS to csv-file.
    FGS_fileName = "FGS_$N_rays.N_subs_$N_subs.Ndim_$Ndim.Kappa_$kappa.Sigmas_$sigma_s.csv"
    CSV.write(FGS_fileName, Tables.table(FGS), header=false)

    # Gas to gas exchange factor matrices.
    # write FSG to csv-file.
    FGG_fileName = "FGG_$N_rays.N_subs_$N_subs.Ndim_$Ndim.Kappa_$kappa.Sigmas_$sigma_s.csv"
    CSV.write(FGG_fileName, Tables.table(FGG), header=false)
    
end