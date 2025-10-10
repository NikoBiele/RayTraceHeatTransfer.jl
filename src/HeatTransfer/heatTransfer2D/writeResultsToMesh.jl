"""
Write results from exchange factor solution back to mesh structure.
Handles grey, spectral_uniform, and spectral_variable modes correctly.
"""

# Grey mode - write scalar results (your existing function with minor fixes)
function write_results_to_mesh_grey!(mesh::RayTracingMeshOptim, T::Vector{P}, j::Vector{P}, 
                                     Abs::Vector{P}, r::Vector{P}, N_surfs::Int) where P
    surf_count = 0
    vol_count = 0
    
    for (i, coarse_face) in enumerate(mesh.coarse_mesh)
        for (k, fine_face) in enumerate(mesh.fine_mesh[i])
            vol_count += 1
            vol_idx = N_surfs + vol_count
            
            # Volume properties - grey mode (scalar)
            e_vol = max(j[vol_idx] - r[vol_idx], 0.0)
            fine_face.T_g = T[vol_idx]
            fine_face.j_g = j[vol_idx]
            fine_face.g_a_g = Abs[vol_idx]
            fine_face.e_g = e_vol
            fine_face.r_g = r[vol_idx]
            fine_face.g_g = Abs[vol_idx] + r[vol_idx]
            fine_face.q_g = e_vol - Abs[vol_idx]

            # Surface properties - grey mode (scalar)
            for (wall_idx, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surf_count += 1
                    e_surf = max(j[surf_count] - r[surf_count], 0.0)
                    
                    fine_face.T_w[wall_idx] = T[surf_count]
                    fine_face.j_w[wall_idx] = j[surf_count]
                    fine_face.g_a_w[wall_idx] = Abs[surf_count]
                    fine_face.e_w[wall_idx] = e_surf
                    fine_face.r_w[wall_idx] = r[surf_count]
                    fine_face.g_w[wall_idx] = Abs[surf_count] + r[surf_count]
                    fine_face.q_w[wall_idx] = e_surf - Abs[surf_count]
                end
            end
        end
    end
    
    println("Grey results written: $vol_count volumes, $surf_count surfaces")
end

# Spectral mode - write results for a specific spectral bin
function write_results_to_mesh_spectral_bin!(mesh::RayTracingMeshOptim, T::Vector{P}, j::Vector{P}, 
                                             Abs::Vector{P}, r::Vector{P}, N_surfs::Int,
                                             spectral_bin::Int) where P
    surf_count = 0
    vol_count = 0
    
    for (i, coarse_face) in enumerate(mesh.coarse_mesh)
        for (k, fine_face) in enumerate(mesh.fine_mesh[i])
            vol_count += 1
            vol_idx = N_surfs + vol_count
            
            # Volume properties - spectral mode (vector)
            e_vol = max(j[vol_idx] - r[vol_idx], 0.0)
            
            # Check if this is first initialization or update
            if isa(fine_face.j_g, Vector)
                # Spectral vector already exists - update specific bin
                fine_face.j_g[spectral_bin] = j[vol_idx]
                fine_face.g_a_g[spectral_bin] = Abs[vol_idx]
                fine_face.e_g[spectral_bin] = e_vol
                fine_face.r_g[spectral_bin] = r[vol_idx]
                fine_face.g_g[spectral_bin] = Abs[vol_idx] + r[vol_idx]
            else
                error("Mesh volume properties are not spectral vectors but spectral_bin=$spectral_bin was specified")
            end
            
            # Temperature and heat source are scalar (updated after all bins)
            fine_face.T_g = T[vol_idx]
            fine_face.q_g = e_vol - Abs[vol_idx]

            # Surface properties - spectral mode (vector per wall)
            for (wall_idx, is_solid) in enumerate(fine_face.solidWalls)
                if is_solid
                    surf_count += 1
                    e_surf = max(j[surf_count] - r[surf_count], 0.0)
                    
                    # Check if this is first initialization or update
                    if isa(fine_face.j_w[wall_idx], Vector)
                        # Spectral vector already exists - update specific bin
                        fine_face.j_w[wall_idx][spectral_bin] = j[surf_count]
                        fine_face.g_a_w[wall_idx][spectral_bin] = Abs[surf_count]
                        fine_face.e_w[wall_idx][spectral_bin] = e_surf
                        fine_face.r_w[wall_idx][spectral_bin] = r[surf_count]
                        fine_face.g_w[wall_idx][spectral_bin] = Abs[surf_count] + r[surf_count]
                    else
                        error("Mesh wall properties are not spectral vectors but spectral_bin=$spectral_bin was specified")
                    end
                    
                    # Temperature and heat source are scalar (updated after all bins)
                    fine_face.T_w[wall_idx] = T[surf_count]
                    fine_face.q_w[wall_idx] = e_surf - Abs[surf_count]
                end
            end
        end
    end
    
    println("Spectral bin $spectral_bin results written: $vol_count volumes, $surf_count surfaces")
end

# Main interface - dispatches based on spectral mode
function write_results_to_mesh!(mesh::RayTracingMeshOptim, T::Vector{P}, j::Vector{P}, 
                               Abs::Vector{P}, r::Vector{P}, N_surfs::Int;
                               spectral_bin::Union{Nothing,Int}=nothing) where P
    
    if mesh.spectral_mode == :grey
        # Grey mode - write scalars
        if spectral_bin !== nothing && spectral_bin != 1
            @warn "spectral_bin=$spectral_bin specified for grey mesh, ignoring"
        end
        write_results_to_mesh_grey!(mesh, T, j, Abs, r, N_surfs)
        
    elseif mesh.spectral_mode == :spectral_uniform || mesh.spectral_mode == :spectral_variable
        # Spectral modes - write to vectors
        if spectral_bin === nothing
            error("spectral_bin must be specified for spectral mesh modes")
        end
        write_results_to_mesh_spectral_bin!(mesh, T, j, Abs, r, N_surfs, spectral_bin)
        
    else
        error("Unknown spectral mode: $(mesh.spectral_mode)")
    end
end

# Convenience function to update temperatures and heat sources after all spectral bins are computed
function update_scalar_temperatures_and_heat_sources!(rtm::RayTracingMeshOptim, wavelength_range::Tuple{Int,Int}=(-7,-3);
                                                      temperatures_in=nothing)
    G = eltype(rtm.fine_mesh[1][1].T_g)

    if rtm.spectral_mode != :spectral_variable
        # Grey or uniform spectral - simple summation
        for ((coarse_idx, fine_idx, wall_index), surface_idx) in rtm.surface_mapping
            sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
            if sub_face.T_in_w[wall_index] < -0.1
                # Radiative equilibrium - compute temperature from spectral emission
                if isa(sub_face.epsilon[wall_index], Vector)
                    epsilon_avg = sum(sub_face.epsilon[wall_index]) / length(sub_face.epsilon[wall_index])
                else
                    epsilon_avg = sub_face.epsilon[wall_index]
                end
                total_emission = sum(sub_face.e_w[wall_index])
                sub_face.T_w[wall_index] = (total_emission / (epsilon_avg * STEFAN_BOLTZMANN * sub_face.area[wall_index]))^(1/4)
            else
                # Prescribed temperature - compute heat source
                sub_face.T_w[wall_index] = sub_face.T_in_w[wall_index]
                sub_face.q_w[wall_index] = sum(sub_face.e_w[wall_index]) - sum(sub_face.g_a_w[wall_index])
            end
        end
        
        for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
            sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
            if sub_face.T_in_g < -0.1
                # Radiative equilibrium - compute temperature from spectral emission
                if isa(sub_face.kappa_g, Vector)
                    kappa_avg = sum(sub_face.kappa_g) / length(sub_face.kappa_g)
                else
                    kappa_avg = sub_face.kappa_g
                end
                total_emission = sum(sub_face.e_g)
                sub_face.T_g = (total_emission / (4 * kappa_avg * STEFAN_BOLTZMANN * sub_face.volume))^(1/4)
            else
                # Prescribed temperature - compute heat source
                sub_face.T_g = sub_face.T_in_g
                sub_face.q_g = sum(sub_face.e_g) - sum(sub_face.g_a_g)
            end
        end
        
    else
        # Spectral variable - use Newton-Raphson for temperature
        if temperatures_in === nothing
            # Get maximum initial temperature
            temperatures_init = zeros(G, length(rtm.surface_mapping) + length(rtm.volume_mapping))
            for ((coarse_idx, fine_idx, wall_index), surface_idx) in rtm.surface_mapping
                sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
                if sub_face.T_in_w[wall_index] > -0.1
                    temperatures_init[surface_idx] = sub_face.T_in_w[wall_index]
                    sub_face.T_w[wall_index] = sub_face.T_in_w[wall_index]
                    # Compute heat source for prescribed temperature
                    sub_face.q_w[wall_index] = sum(sub_face.e_w[wall_index]) - sum(sub_face.g_a_w[wall_index])
                end
            end
            for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
                sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
                if sub_face.T_in_g > -0.1
                    temperatures_init[length(rtm.surface_mapping) + volume_idx] = sub_face.T_in_g
                    sub_face.T_g = sub_face.T_in_g
                    # Compute heat source for prescribed temperature
                    sub_face.q_g = sum(sub_face.e_g) - sum(sub_face.g_a_g)
                end
            end
            maximum_temperature = maximum(temperatures_init)
            
            # Solve for temperatures using Newton-Raphson
            wavelength_bands = G.(wavelength_band_splits(rtm, wavelength_range))
            for ((coarse_idx, fine_idx, wall_index), surface_idx) in rtm.surface_mapping
                sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
                if sub_face.T_in_w[wall_index] < -0.1
                    sub_face.T_w[wall_index] = solve_temperature_newton_raphson(
                        rtm, sub_face.area[wall_index], sub_face.e_w[wall_index],
                        sub_face.epsilon[wall_index], wavelength_bands; 
                        initial_temp=maximum_temperature, max_iter=10_000, tolerance=1000*eps(G)
                    )
                    # Heat source is computed from energy balance
                    sub_face.q_w[wall_index] = sum(sub_face.e_w[wall_index]) - sum(sub_face.g_a_w[wall_index])
                end
            end
            for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
                sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
                if sub_face.T_in_g < -0.1
                    sub_face.T_g = solve_temperature_newton_raphson(
                        rtm, 4*sub_face.volume, sub_face.e_g,
                        sub_face.kappa_g, wavelength_bands; 
                        initial_temp=maximum_temperature, max_iter=10_000, tolerance=1000*eps(G)
                    )
                    # Heat source is computed from energy balance
                    sub_face.q_g = sum(sub_face.e_g) - sum(sub_face.g_a_g)
                end
            end
        end
    end
end