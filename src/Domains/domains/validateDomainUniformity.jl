function validateSpectralUniformity!(rtm::RayTracingDomain2D; atol=1e-10, verbose=true)
    """
    Checks if spectral properties are uniform across all faces and spectral bins.
    """
    first_epsilon = nothing
    first_kappa_g = nothing
    first_sigma_s_g = nothing
    
    for ((coarse_face, fine_face, wall_index), surface_index) in rtm.surface_mapping
        face = rtm.fine_mesh[coarse_face][fine_face]
        first_epsilon = nothing
        for bin in 1:rtm.n_spectral_bins
            if bin == 1
                first_epsilon = face.epsilon[wall_index][bin]
            else
                epsilon = face.epsilon[wall_index][bin]
                if abs(epsilon - first_epsilon) > atol
                    verbose && println("Spectral variation detected across walls, using spectral solver")
                    return false
                end
            end
        end        
    end
    verbose && println("No spectral variation detected separately across walls")
    
    for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping        
        face = rtm.fine_mesh[coarse_face][fine_face]
        first_kappa_g = nothing
        first_sigma_s_g = nothing
        all(isfinite, face.kappa_g) && all(isfinite, face.sigma_s_g) ||
            error("Non-finite spectral properties on face ($coarse_face, $fine_face)")
            
        for bin in 1:rtm.n_spectral_bins
            if bin == 1
                first_kappa_g = face.kappa_g[bin]
                first_sigma_s_g = face.sigma_s_g[bin]
            else
                kappa_g = face.kappa_g[bin]
                sigma_s_g = face.sigma_s_g[bin]
                if abs(kappa_g - first_kappa_g) > atol || abs(sigma_s_g - first_sigma_s_g) > atol
                    verbose && println("Spectral variation detected across volumes, using spectral solver")
                    return false
                end
            end
        end
    end
    verbose && println("No spectral variation detected separately across volumes")
    
    use_efficient = isapprox(first_epsilon, first_kappa_g/(first_kappa_g+first_sigma_s_g), rtol=1e-10) && 
                    abs(first_epsilon - 1.0) < 1e-10 ? true : false
    (verbose && use_efficient) && println("No spectral variation detected across mesh, using efficient grey solver")
    (verbose && !use_efficient) && println("Spectral variation detected across mesh, using full spectral solver")

    return use_efficient
end

function validateExtinctionUniformity!(rtm::RayTracingDomain2D; atol=1e-5, verbose=true)
    """
    Checks if extinction properties are uniform across all faces and spectral bins.
    """

    # bin the spectral bins into uniform extinction blocks
    uniform_across_bin = []
    for bin = 1:rtm.n_spectral_bins
        bin_break = false
        first_beta_bin = nothing
        for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping
            face = rtm.fine_mesh[coarse_face][fine_face]
            if first_beta_bin === nothing
                first_beta_bin = face.kappa_g[bin] + face.sigma_s_g[bin]
            else
                if abs(first_beta_bin - (face.kappa_g[bin] + face.sigma_s_g[bin])) > atol
                    push!(uniform_across_bin, -1.0)
                    bin_break = true
                    break
                end
            end
        end
        if !bin_break
            push!(uniform_across_bin, first_beta_bin)
        end
    end
    
    return uniform_across_bin
end

# Domain uniformity validation functions for 3D
# These check if properties are uniform across the domain to enable optimizations

function inheritSpectralUniformity3D!(rtm::RayTracingDomain2D) #; atol=1e-5)
    """
    Checks if spectral properties are uniform across all faces and spectral bins.
    """
    return rtm.spectral_mode
    
end

function inheritExtinctionUniformity3D!(rtm::RayTracingDomain2D) # rtm::RayTracingDomain3D; atol=1e-5)
    """
    Checks if extinction properties are uniform across all volumes and spectral bins.
    This enables faster ray tracing using a single beta value.
    """
    return rtm.uniform_across_bin

end