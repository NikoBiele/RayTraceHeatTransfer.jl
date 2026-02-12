# Main interface - dispatches based on spectral mode
function writeResultsToDomain!(domain::ViewFactorDomain3D, j::Vector{P}, 
                                    Abs::Vector{P}, r::Vector{P};
                                    T::Union{Nothing,Vector{P}}=nothing, spectral_bin::Union{Nothing,Int}=nothing) where P
    
    if domain.spectral_mode == :grey
        # Grey mode
        if spectral_bin !== nothing && spectral_bin != 1
            @warn "spectral_bin=$spectral_bin specified for grey domain, ignoring"
        end
        writeResultsToDomainGrey3D!(domain, T, j, Abs, r)
        
    elseif domain.spectral_mode == :spectral_uniform || domain.spectral_mode == :spectral_variable
        # Spectral mode
        if spectral_bin === nothing
            error("spectral_bin must be specified for spectral domain")
        end
        writeResultsToDomainSpectralBin3D!(domain, j, Abs, r, spectral_bin)
        
    else
        error("Unknown spectral mode: $(domain.spectral_mode)")
    end
end

# Main interface - dispatches based on spectral mode
function writeResultsToDomain!(domain::RayTracingDomain2D, j::Vector{P}, 
                               Abs::Vector{P}, r::Vector{P};
                               T::Union{Nothing,Vector{P}}=nothing, spectral_bin::Union{Nothing,Int}=nothing) where P
    
    if domain.spectral_mode == :grey
        # Grey mode - write scalars
        if spectral_bin !== nothing && spectral_bin != 1
            @warn "spectral_bin=$spectral_bin specified for grey domain, ignoring"
        end
        writeResultsToDomainGrey!(domain, T, j, Abs, r)
        
    elseif domain.spectral_mode == :spectral_uniform || domain.spectral_mode == :spectral_variable
        # Spectral modes - write to vectors
        if spectral_bin === nothing
            error("spectral_bin must be specified for spectral domain modes")
        end
        writeResultsToDomainSpectralBin!(domain, j, Abs, r, spectral_bin)
        
    else
        error("Unknown spectral mode: $(domain.spectral_mode)")
    end
end

# Grey mode - write scalar results
function writeResultsToDomainGrey3D!(domain::ViewFactorDomain3D, T::Vector{P}, j::Vector{P}, 
                                          Abs::Vector{P}, r::Vector{P}) where P
    surf_count = 0
    
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            
            e_surf = max(j[surf_count] - r[surf_count], 0.0)
            
            # Grey mode - scalar values
            subface.j_w = j[surf_count]
            subface.g_a_w = Abs[surf_count]
            subface.e_w = e_surf
            subface.r_w = r[surf_count]
            subface.g_w = Abs[surf_count] + r[surf_count]
            subface.i_w = j[surf_count] / (π * subface.area)
            if subface.T_in_w < -0.1
                # Prescribed flux
                subface.q_w = subface.q_in_w
                subface.T_w = T[surf_count]
            else
                subface.T_w = subface.T_in_w
                subface.q_w = e_surf - Abs[surf_count]
            end
        end
    end
    
    println("Grey results written: $surf_count surfaces")
end

# Spectral mode - write results for a specific spectral bin
function writeResultsToDomainSpectralBin3D!(domain::ViewFactorDomain3D, j::Vector{P}, 
                                                  Abs::Vector{P}, r::Vector{P}, spectral_bin::Int) where P
    surf_count = 0
    
    for superface in domain.facesMesh
        for subface in superface.subFaces
            surf_count += 1
            
            e_surf = max(j[surf_count] - r[surf_count], 0.0)
            
            # Check if spectral vectors already initialized
            if isa(subface.j_w, Vector)
                # Update specific bin
                subface.j_w[spectral_bin] = j[surf_count]
                subface.g_a_w[spectral_bin] = Abs[surf_count]
                subface.e_w[spectral_bin] = e_surf
                subface.r_w[spectral_bin] = r[surf_count]
                subface.g_w[spectral_bin] = Abs[surf_count] + r[surf_count]
                subface.i_w[spectral_bin] = j[surf_count] / (π * subface.area)
            else
                error("Domain surface properties are not spectral vectors but spectral_bin=$spectral_bin was specified")
            end
            
        end
    end
    
    println("Spectral bin $spectral_bin results written: $surf_count surfaces")
end

# Grey mode - write scalar results (your existing function with minor fixes)
function writeResultsToDomainGrey!(domain::RayTracingDomain2D, T::Vector{P}, j::Vector{P}, 
                                     Abs::Vector{P}, r::Vector{P}) where P

    for ((coarse_face, fine_face_idx, wall_idx), surface_idx) in domain.surface_mapping
        fine_face = domain.fine_mesh[coarse_face][fine_face_idx]
        e_surf = max(j[surface_idx] - r[surface_idx], 0.0)

        fine_face.T_w[wall_idx] = T[surface_idx]
        fine_face.j_w[wall_idx] = j[surface_idx]
        fine_face.g_a_w[wall_idx] = Abs[surface_idx]
        fine_face.e_w[wall_idx] = e_surf
        fine_face.r_w[wall_idx] = r[surface_idx]
        fine_face.g_w[wall_idx] = Abs[surface_idx] + r[surface_idx]
        fine_face.q_w[wall_idx] = e_surf - Abs[surface_idx]
    end
    N_surfs = length(domain.surface_mapping)
    for ((coarse_face, fine_face_idx), volume_idx) in domain.volume_mapping
        fine_face = domain.fine_mesh[coarse_face][fine_face_idx]
        vol_idx = N_surfs + volume_idx    
        e_vol = max(j[vol_idx] - r[vol_idx], 0.0)
        fine_face.T_g = T[vol_idx]
        fine_face.j_g = j[vol_idx]
        fine_face.g_a_g = Abs[vol_idx]
        fine_face.e_g = e_vol
        fine_face.r_g = r[vol_idx]
        fine_face.g_g = Abs[vol_idx] + r[vol_idx]
        fine_face.q_g = e_vol - Abs[vol_idx]
    end
end

# Spectral mode - write results for a specific spectral bin
function writeResultsToDomainSpectralBin!(domain::RayTracingDomain2D, j::Vector{P}, 
                                             Abs::Vector{P}, r::Vector{P},
                                             spectral_bin::Int) where P

    for ((coarse_face, fine_face_idx, wall_idx), surface_idx) in domain.surface_mapping
        fine_face = domain.fine_mesh[coarse_face][fine_face_idx]
        e_surf = max(j[surface_idx] - r[surface_idx], 0.0)

        fine_face.j_w[wall_idx][spectral_bin] = j[surface_idx]
        fine_face.g_a_w[wall_idx][spectral_bin] = Abs[surface_idx]
        fine_face.e_w[wall_idx][spectral_bin] = e_surf
        fine_face.r_w[wall_idx][spectral_bin] = r[surface_idx]
        fine_face.g_w[wall_idx][spectral_bin] = Abs[surface_idx] + r[surface_idx]
    end
    N_surfs = length(domain.surface_mapping)
    for ((coarse_face, fine_face_idx), volume_idx) in domain.volume_mapping
        fine_face = domain.fine_mesh[coarse_face][fine_face_idx]
        vol_idx = N_surfs + volume_idx    
        e_vol = max(j[vol_idx] - r[vol_idx], 0.0)
        fine_face.j_g[spectral_bin] = j[vol_idx]
        fine_face.g_a_g[spectral_bin] = Abs[vol_idx]
        fine_face.e_g[spectral_bin] = e_vol
        fine_face.r_g[spectral_bin] = r[vol_idx]
        fine_face.g_g[spectral_bin] = Abs[vol_idx] + r[vol_idx]
    end
end