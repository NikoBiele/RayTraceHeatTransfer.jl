"""
    GasProperties

This struct holds information about the radiative properties of
the participating medium.
"""
struct GasProperties
    sigma_s::Float64 # scattering coefficient
    kappa::Float64 # absorption coefficient
    beta::Float64 # extinction coefficient
    omega::Float64 # scattering albedo
    function GasProperties(sigma_s::Float64,kappa::Float64)
        beta = sigma_s+kappa
        omega = sigma_s/beta 
        return new(sigma_s,kappa,beta,omega)
    end
end