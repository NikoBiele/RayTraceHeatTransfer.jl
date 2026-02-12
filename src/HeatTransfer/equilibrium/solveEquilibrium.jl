function solveEquilibrium!(domain::Union{RayTracingDomain2D, ViewFactorDomain3D}, F_matrices::Union{Matrix{G}, Vector{Matrix{G}}};
                            max_iterations::Int=500, spectral_coupling_tolerance::G=sqrt(eps(G))) where G
    if domain isa RayTracingDomain2D
        if domain.spectral_mode == :grey
            # Grey mode
            equilibriumGrey2D!(domain, F_matrices)
        elseif domain.spectral_mode == :spectral_uniform || domain.spectral_mode == :spectral_variable
            # Spectral mode
            equilibriumSpectral2D!(domain, F_matrices; max_iterations=max_iterations)
        else
            error("Unknown spectral mode: $(domain.spectral_mode)")
        end
    elseif domain isa ViewFactorDomain3D
        if domain.spectral_mode == :grey
            # Grey mode
            equilibriumSurfacesGrey3D!(domain, F_matrices)
        elseif domain.spectral_mode == :spectral_uniform || domain.spectral_mode == :spectral_variable
            # Spectral mode
            equilibriumSurfacesSpectral3D!(domain, F_matrices; max_iterations=max_iterations)
        else
            error("Unknown spectral mode: $(domain.spectral_mode)")
        end
    else
        error("Unknown domain type: $(typeof(domain))")
    end
end