function validateSpectralUniformity!(rtm::RayTracingDomain2D; atol=1e-5)
    """
    Checks if spectral properties are uniform across all faces and spectral bins.
    """
    
    first_epsilon = nothing
    for ((coarse_face, fine_face, wall_index), surface_index) in rtm.surface_mapping
        face = rtm.fine_mesh[coarse_face][fine_face]
        for bin in 1:rtm.n_spectral_bins
            if first_epsilon === nothing
                first_epsilon = face.epsilon[wall_index][bin]
            else
                epsilon = face.epsilon[wall_index][bin]
                if abs(epsilon - first_epsilon) > atol
                    println("Spectral variation detected across mesh, using spectral solver")
                    return false
                end
            end
        end        
    end
    println("No spectral variation detected across walls")
    
    first_kappa_g = nothing
    first_sigma_s_g = nothing
    for ((coarse_face, fine_face), volume_index) in rtm.volume_mapping        
        face = rtm.fine_mesh[coarse_face][fine_face]
        for bin in 1:rtm.n_spectral_bins
            if first_kappa_g === nothing
                first_kappa_g = face.kappa_g[bin]
                first_sigma_s_g = face.sigma_s_g[bin]
            else
                kappa_g = face.kappa_g[bin]
                sigma_s_g = face.sigma_s_g[bin]
                if abs(kappa_g - first_kappa_g) > atol || abs(sigma_s_g - first_sigma_s_g) > atol
                    println("Spectral variation detected across mesh, using spectral solver")
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
    
    # first_epsilon = nothing
    # for ((coarse_vol, fine_vol, depth_vol, face_index), surface_index) in rtm.surface_mapping
    #     volume = rtm.fine_mesh[coarse_vol][fine_vol][depth_vol]
    #     face = volume.faces[face_index]
        
    #     for bin in 1:rtm.n_spectral_bins
    #         if first_epsilon === nothing
    #             if isa(face.epsilon, Vector)
    #                 first_epsilon = face.epsilon[bin]
    #             else
    #                 first_epsilon = face.epsilon
    #             end
    #         else
    #             epsilon = isa(face.epsilon, Vector) ? face.epsilon[bin] : face.epsilon
    #             if abs(epsilon - first_epsilon) > atol
    #                 println("Spectral variation detected across mesh, using spectral solver")
    #                 return false
    #             end
    #         end
    #     end        
    # end
    # println("No spectral variation detected across faces")
    
    # first_kappa_g = nothing
    # first_sigma_s_g = nothing
    # for ((coarse_vol, fine_vol, depth_vol), volume_index) in rtm.volume_mapping        
    #     volume = rtm.fine_mesh[coarse_vol][fine_vol][depth_vol]
        
    #     for bin in 1:rtm.n_spectral_bins
    #         if first_kappa_g === nothing
    #             if isa(volume.kappa_g, Vector)
    #                 first_kappa_g = volume.kappa_g[bin]
    #                 first_sigma_s_g = volume.sigma_s_g[bin]
    #             else
    #                 first_kappa_g = volume.kappa_g
    #                 first_sigma_s_g = volume.sigma_s_g
    #             end
    #         else
    #             kappa_g = isa(volume.kappa_g, Vector) ? volume.kappa_g[bin] : volume.kappa_g
    #             sigma_s_g = isa(volume.sigma_s_g, Vector) ? volume.sigma_s_g[bin] : volume.sigma_s_g
                
    #             if abs(kappa_g - first_kappa_g) > atol || abs(sigma_s_g - first_sigma_s_g) > atol
    #                 println("Spectral variation detected across mesh, using spectral solver")
    #                 return false
    #             end
    #         end
    #     end
    # end
    # println("No spectral variation detected across volumes")
    
    # println("No spectral variation detected across mesh, using efficient grey solver")
    # return true
end

function inheritExtinctionUniformity3D!(rtm::RayTracingDomain2D) # rtm::RayTracingDomain3D; atol=1e-5)
    """
    Checks if extinction properties are uniform across all volumes and spectral bins.
    This enables faster ray tracing using a single beta value.
    """
    return rtm.uniform_extinction

    # first_beta = nothing
    
    # for ((coarse_vol, fine_vol, depth_vol), volume_index) in rtm.volume_mapping
    #     volume = rtm.fine_mesh[coarse_vol][fine_vol][depth_vol]

    #     for bin = 1:rtm.n_spectral_bins
    #         if first_beta === nothing
    #             if isa(volume.kappa_g, Vector)
    #                 first_beta = volume.kappa_g[bin] + volume.sigma_s_g[bin]
    #             else
    #                 first_beta = volume.kappa_g + volume.sigma_s_g
    #             end
    #         else
    #             beta = if isa(volume.kappa_g, Vector)
    #                 volume.kappa_g[bin] + volume.sigma_s_g[bin]
    #             else
    #                 volume.kappa_g + volume.sigma_s_g
    #             end
                
    #             if abs(first_beta - beta) > atol
    #                 println("Extinction variation detected across the spectrum, ray tracing each spectral bin separately")
    #                 return false
    #             end
    #         end
    #     end
    # end

    # println("No extinction variation detected across the spectrum, ray tracing grey domain only")
    # return true
end