function validateSpectralUniformity!(rtm::RayTracingDomain2D; atol=1e-5)
    """
    Checks if spectral properties are uniform across all faces and spectral bins.
    """
    
    for ((coarse_face, fine_face, wall_index), surface_index) in rtm.surface_mapping
        face = rtm.fine_mesh[coarse_face][fine_face]
        first_epsilon = nothing
        for bin in 1:rtm.n_spectral_bins
            if bin == 1
                first_epsilon = face.epsilon[wall_index][bin]
            else
                epsilon = face.epsilon[wall_index][bin]
                if abs(epsilon - first_epsilon) > atol
                    println("Spectral variation detected across walls, using spectral solver")
                    return false
                end
            end
        end        
    end
    println("No spectral variation detected across walls")
    
    for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping        
        face = rtm.fine_mesh[coarse_face][fine_face]
        first_kappa_g = nothing
        first_sigma_s_g = nothing
        for bin in 1:rtm.n_spectral_bins
            if bin == 1
                first_kappa_g = face.kappa_g[bin]
                first_sigma_s_g = face.sigma_s_g[bin]
            else
                kappa_g = face.kappa_g[bin]
                sigma_s_g = face.sigma_s_g[bin]
                if abs(kappa_g - first_kappa_g) > atol || abs(sigma_s_g - first_sigma_s_g) > atol
                    println("Spectral variation detected across volumes, using spectral solver")
                    return false
                end
            end
        end
    end
    println("No spectral variation detected across volumes")
    
    println("No spectral variation detected across mesh, using efficient grey solver")
    return true
end

function validateExtinctionUniformity!(rtm::RayTracingDomain2D; atol=1e-5)
    """
    Checks if extinction properties are uniform across all faces and spectral bins.
    """
    first_beta = nothing
    for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping
        face = rtm.fine_mesh[coarse_face][fine_face]
        for bin = 1:rtm.n_spectral_bins
            if first_beta === nothing
                first_beta = face.kappa_g[bin] + face.sigma_s_g[bin]
            else
                if abs(first_beta - (face.kappa_g[bin] + face.sigma_s_g[bin])) > atol
                    println("Extinction variation detected across the spectrum, ray tracing each spectral bin separately")
                    return false
                end
            end
        end
    end
    
    println("No extinction variation detected across the spectrum, ray tracing grey domain only")
    return true
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
    return rtm.uniform_extinction

end