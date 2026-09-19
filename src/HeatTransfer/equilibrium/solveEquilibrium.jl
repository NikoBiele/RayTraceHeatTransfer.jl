function solveEquilibrium!(domain::Union{RayTracingDomain2D, SurfaceDomain3D}, F_matrices::Union{AbstractMatrix, Vector{AbstractMatrix}};
                            max_iters::Int=1000, convergence_tol::G=1e-12, verbose::Bool=true) where G
    if domain isa RayTracingDomain2D
        if domain.directional_model !== nothing
            if domain.spectral_mode == :grey
                return equilibriumGreyDirectional2D!(domain, F_matrices; verbose=verbose)
            elseif domain.spectral_mode == :spectral_variable
                return equilibriumSpectralDirectional2D!(domain, F_matrices; max_iters=max_iters,
                                                         convergence_tol=convergence_tol, verbose=verbose)
            else
                error("the domain is :spectral_uniform (no reflection or scattering), so a directional_model " *
                      "has nothing to redistribute; set domain.directional_model = nothing")
            end
        end
        if domain.spectral_mode == :grey
            # Grey mode
            equilibriumGrey2D!(domain, F_matrices; verbose=verbose)
        elseif domain.spectral_mode == :spectral_uniform || domain.spectral_mode == :spectral_variable
            # Spectral mode
            equilibriumSpectral2D!(domain, F_matrices; max_iters=max_iters, convergence_tol=convergence_tol, verbose=verbose)
        else
            error("Unknown spectral mode: $(domain.spectral_mode)")
        end
    elseif domain isa SurfaceDomain3D
        if domain.spectral_mode == :grey
            # Grey mode
            equilibriumSurfacesGrey3D!(domain, F_matrices; verbose=verbose)
        elseif domain.spectral_mode == :spectral_uniform || domain.spectral_mode == :spectral_variable
            # Spectral mode
            equilibriumSurfacesSpectral3D!(domain, F_matrices; max_iters=max_iters, convergence_tol=convergence_tol, verbose=verbose)
        else
            error("Unknown spectral mode: $(domain.spectral_mode)")
        end
    else
        error("Unknown domain type: $(typeof(domain))")
    end
end