function wavelength_band_splits(rtm::RayTracingMeshOptim, wavelength_range::Tuple{Int,Int}=(-1,-1))
    # Determine spectral configuration
    n_bins = rtm.n_spectral_bins
    
    if wavelength_range[1] >= 0 || wavelength_range[2] >= 0
        @warn "Wavelength range is invalid or not specified. Defaulting to 1e-12 to 1e3"
        wavelength_bands = setupWavelengthBands()
    else
        if length(wavelength_range) != 2
            error("wavelength_range must be a tuple of length 2")
        end
        if wavelength_range[1] > wavelength_range[2]
            error("wavelength_range[1] must be less than wavelength_range[2]")
        end
        wavelength_bands = setupWavelengthBands(wavelength_range[1], wavelength_range[2], n_bins)
    end
    
    return wavelength_bands
end