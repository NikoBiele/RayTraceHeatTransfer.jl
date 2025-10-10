function setupWavelengthBands(log10_min_split::Int=-7, log10_max_split::Int=-3, n_bins::Int=10)
   
    # Default wavelength splits - you can customize this
    wavelength_bands = 10.0 .^range(log10_min_split, log10_max_split, length=n_bins-2)
    wavelength_bands = [1e-8; wavelength_bands; 1.0] # from 'zero' to 'infinity'

    return wavelength_bands
end