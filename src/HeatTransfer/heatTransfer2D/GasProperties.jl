struct GasProperties{T}
    kappa::T # absorption coefficient
    sigma_s::T # scattering coefficient
    beta::T # extinction coefficient
    omega::T # single scattering albedo
end
function GasProperties(kappa::T,sigma_s::T) where T
    beta = sigma_s+kappa
    omega = sigma_s/beta 
    return GasProperties(kappa,sigma_s,beta,omega)
end