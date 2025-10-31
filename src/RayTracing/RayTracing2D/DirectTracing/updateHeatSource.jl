function update_spectral_results!(rtm::RayTracingDomain2D, absorbed_count::Vector{Vector{P}},
                            gas_emitted_count::Vector{Vector{P}}, wall_emitted_count::Vector{Vector{Vector{P}}},
                            reflected_count::Vector{Vector{Vector{P}}}, scattered_count::Vector{Vector{P}},
                            wall_absorbed_count::Vector{Vector{Vector{P}}}, total_energy::G, num_rays::P,
                            spectral_bin::P=1) where {G, P<:Integer}

    energy_per_ray = G(total_energy / num_rays)

    println("Updating spectral results for spectral bin $spectral_bin")
    println("Energy per ray: $energy_per_ray")
    
    for ((coarse_index, fine_index, wall_index), surface_index) in rtm.surface_mapping
        sub_face = rtm.fine_mesh[coarse_index][fine_index]
        wall_absorbed = G(wall_absorbed_count[coarse_index][fine_index][wall_index] * energy_per_ray)
        wall_emitted = G(wall_emitted_count[coarse_index][fine_index][wall_index] * energy_per_ray)
        wall_reflected = G(reflected_count[coarse_index][fine_index][wall_index] * energy_per_ray)
        
        # Update wall power variables - handle spectral vs grey
        if isa(sub_face.g_a_w[wall_index], Vector)
            # Spectral case
            sub_face.g_a_w[wall_index][spectral_bin] = wall_absorbed
            sub_face.e_w[wall_index][spectral_bin] = wall_emitted
            sub_face.r_w[wall_index][spectral_bin] = wall_reflected
            sub_face.j_w[wall_index][spectral_bin] = wall_emitted + wall_reflected
            sub_face.g_w[wall_index][spectral_bin] = wall_absorbed + wall_reflected
        else
            # Grey case
            if spectral_bin == 1
                sub_face.g_a_w[wall_index] = wall_absorbed
                sub_face.e_w[wall_index] = wall_emitted
                sub_face.r_w[wall_index] = wall_reflected
                sub_face.j_w[wall_index] = wall_emitted + wall_reflected
                sub_face.g_w[wall_index] = wall_absorbed + wall_reflected
            end
        end
    end

    for ((coarse_index, fine_index), volume_index) in rtm.volume_mapping
        sub_face = rtm.fine_mesh[coarse_index][fine_index]

        # Gas properties for this spectral bin
        gas_absorbed = G(absorbed_count[coarse_index][fine_index] * energy_per_ray)
        gas_emitted = G(gas_emitted_count[coarse_index][fine_index] * energy_per_ray)
        gas_scattered = G(scattered_count[coarse_index][fine_index] * energy_per_ray)
        
        # Update gas properties - handle spectral vs grey
        if isa(sub_face.g_a_g, Vector)
            # Spectral case - update specific bin
            sub_face.g_a_g[spectral_bin] = gas_absorbed
            sub_face.e_g[spectral_bin] = gas_emitted
            sub_face.r_g[spectral_bin] = gas_scattered
            sub_face.j_g[spectral_bin] = gas_emitted + gas_scattered
            sub_face.g_g[spectral_bin] = gas_absorbed + gas_scattered
        else
            # Grey case - direct assignment (should only happen for bin 1)
            if spectral_bin == 1
                sub_face.g_a_g = gas_absorbed
                sub_face.e_g = gas_emitted
                sub_face.r_g = gas_scattered
                sub_face.j_g = gas_emitted + gas_scattered
                sub_face.g_g = gas_absorbed + gas_scattered
            end
        end 
    end
end

function update_scalar_temperatures_and_heat_sources_direct!(rtm::RayTracingDomain2D;
                                                    temperatures_in=nothing)
    G = eltype(rtm.fine_mesh[1][1].T_g)

    if rtm.spectral_mode != :spectral_variable
        # grey or uniform spectral
        temperatures_calc = zeros(G, length(rtm.surface_mapping) + length(rtm.volume_mapping))
        for ((coarse_idx, fine_idx, wall_index), surface_idx) in rtm.surface_mapping
            sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
            if sub_face.T_in_w[wall_index] < -0.1
                epsilon_wall = sum(sub_face.epsilon[wall_index])/length(sub_face.epsilon[wall_index])
                sub_face.T_w[wall_index] = (sum(sub_face.e_w[wall_index]) / (epsilon_wall * STEFAN_BOLTZMANN * sub_face.area[wall_index]))^(1/4)
                temperatures_calc[surface_idx] = sub_face.T_w[wall_index]
            else
                sub_face.T_w[wall_index] = sub_face.T_in_w[wall_index]
                temperatures_calc[surface_idx] = sub_face.T_w[wall_index]
                # update heat source
                sub_face.q_w[wall_index] = sum(sub_face.e_w[wall_index]) - sum(sub_face.g_a_w[wall_index])
            end
        end
        for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
            sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
            if sub_face.T_in_g < -0.1
                local_kappa = sum(sub_face.kappa_g)/length(sub_face.kappa_g)
                sub_face.T_g = (sum(sub_face.e_g) / (4 * local_kappa * STEFAN_BOLTZMANN * sub_face.volume))^(1/4)
                temperatures_calc[length(rtm.surface_mapping) + volume_idx] = sub_face.T_g
            else
                sub_face.T_g = sub_face.T_in_g
                temperatures_calc[length(rtm.surface_mapping) + volume_idx] = sub_face.T_g
                # update heat source
                sub_face.q_g = sum(sub_face.e_g) - sum(sub_face.g_a_g)
            end
        end
    else # spectral variable
        if temperatures_in === nothing
            # first get the maximum initial temperature
            temperatures_init = zeros(G, length(rtm.surface_mapping) + length(rtm.volume_mapping))
            for ((coarse_idx, fine_idx, wall_index), surface_idx) in rtm.surface_mapping
                sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
                if sub_face.T_in_w[wall_index] > -0.1
                    temperatures_init[surface_idx] = sub_face.T_in_w[wall_index]
                    sub_face.T_w[wall_index] = sub_face.T_in_w[wall_index]
                end
            end
            for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
                sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
                if sub_face.T_in_g > -0.1
                    temperatures_init[length(rtm.surface_mapping) + volume_idx] = sub_face.T_in_g
                    sub_face.T_g = sub_face.T_in_g
                end
            end
            maximum_temperature = maximum(temperatures_init)
            # next, find the matching temperatures
            for ((coarse_idx, fine_idx, wall_index), surface_idx) in rtm.surface_mapping
                sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
                if sub_face.T_in_w[wall_index] < -0.1
                    sub_face.T_w[wall_index] = solve_temperature_newton_raphson(rtm, sub_face.area[wall_index], sub_face.e_w[wall_index],
                                            sub_face.epsilon[wall_index]; 
                                            initial_temp=maximum_temperature, max_iter=10_000, tolerance=sqrt(eps(G)))
                end
            end
            for ((coarse_idx, fine_idx), volume_idx) in rtm.volume_mapping
                sub_face = rtm.fine_mesh[coarse_idx][fine_idx]
                if sub_face.T_in_g < -0.1
                    sub_face.T_g = solve_temperature_newton_raphson(rtm, 4*sub_face.volume, sub_face.e_g,
                                            sub_face.kappa_g; 
                                            initial_temp=maximum_temperature, max_iter=10_000, tolerance=sqrt(eps(G)))
                end
            end
        end
    end
end