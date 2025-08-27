function update_heat_source!(rtm::RayTracingMeshOptim, absorbed_count::Vector{Vector{Int}},
                            gas_emitted_count::Vector{Vector{Int}}, wall_emitted_count::Vector{Vector{Vector{Int}}},
                            reflected_count::Vector{Vector{Vector{Int}}}, scattered_count::Vector{Vector{Int}},
                            wall_absorbed_count::Vector{Vector{Vector{Int}}}, total_energy::Float64, num_rays::Int, gas::GasProperties)

    energy_per_ray = total_energy / num_rays

    println("Updating heat source...")
    println("Energy per ray: $energy_per_ray")
    
    total_gas_elements = 0
    nonzero_temps = 0
    equilibrium_elements = 0
    max_temp = 0.0
    total_gas_absorbed_power = 0.0
    total_gas_emitted_power = 0.0
    
    for (coarse_index, coarse_face) in enumerate(rtm.coarse_mesh)
        fine_index = 1

        for sub_face in rtm.fine_mesh[coarse_index]
            total_gas_elements += 1
            
            # Gas properties
            gas_absorbed = absorbed_count[coarse_index][fine_index] * energy_per_ray
            gas_emitted = gas_emitted_count[coarse_index][fine_index] * energy_per_ray
            gas_scattered = scattered_count[coarse_index][fine_index] * energy_per_ray
            
            total_gas_absorbed_power += gas_absorbed
            total_gas_emitted_power += gas_emitted
            
            # Update gas properties from ray tracing results
            sub_face.g_a_g = gas_absorbed  # Absorbed power for gas  
            sub_face.e_g = gas_emitted     # Emitted power for gas (from ray tracing)
            sub_face.r_g = gas_scattered   # Scattered power for gas
            
            # For radiative equilibrium: emitted should equal absorbed
            # The ray tracing gives us the transport, but for equilibrium we need emission = absorption
            sub_face.j_g = sub_face.e_g + sub_face.r_g     # Total outgoing radiant power
            sub_face.g_g = sub_face.g_a_g + sub_face.r_g   # Total incident power (for energy balance)
            
            # Heat source calculation
            if sub_face.T_in_g < 0.0
                # Radiative equilibrium: q_g = 0 (no net energy gain/loss)
                sub_face.q_g = 0.0
            else
                # Boundary element: q_g = absorbed - emitted (net energy balance)
                sub_face.q_g = gas_emitted - gas_absorbed
            end
            
            # Temperature calculation logic based on your convention (check input field)
            if sub_face.T_in_g < 0.0
                equilibrium_elements += 1
                # Radiative equilibrium element (negative T_in indicates unknown temperature)
                # For equilibrium: absorbed = emitted, so calculate T from emitted power
                if gas_emitted > 0.0 && sub_face.volume > 0.0 && gas.kappa > 0.0
                    # T calculated from emitted power: e = 4*kappa*sigma*T^4*volume
                    temp_calc = gas_emitted / (4 * gas.kappa * STEFAN_BOLTZMANN * sub_face.volume)
                    if temp_calc > 0.0
                        sub_face.T_g = temp_calc^(1/4)  # Calculate and store in T_g
                        nonzero_temps += 1
                        max_temp = max(max_temp, sub_face.T_g)
                        
                        # Print first few calculated temperatures
                        if nonzero_temps <= 5
                            println("Calculated T_g for equilibrium element ($coarse_index,$fine_index): $(sub_face.T_g) K, emitted = $gas_emitted W, absorbed = $gas_absorbed W")
                        end
                    else
                        # No emission yet, keep T_g as is or set to small value
                        if nonzero_temps <= 5
                            println("Element ($coarse_index,$fine_index): No emission yet, T_in_g = $(sub_face.T_in_g)")
                        end
                    end
                else
                    if nonzero_temps <= 5
                        println("Element ($coarse_index,$fine_index): Cannot calculate temperature, emitted = $gas_emitted W, T_in_g = $(sub_face.T_in_g)")
                    end
                end
            else
                # Prescribed boundary temperature (non-negative T_in) - copy to T_g
                sub_face.T_g = sub_face.T_in_g
                if nonzero_temps <= 5  # Print first few boundary elements too
                    println("Prescribed temperature element ($coarse_index,$fine_index): T_g = $(sub_face.T_g) K (from T_in_g)")
                end
            end
            
            # Update wall properties
            for wall_index in 1:length(sub_face.solidWalls)
                if sub_face.solidWalls[wall_index]
                    wall_absorbed = wall_absorbed_count[coarse_index][fine_index][wall_index] * energy_per_ray
                    wall_emitted = wall_emitted_count[coarse_index][fine_index][wall_index] * energy_per_ray
                    wall_reflected = reflected_count[coarse_index][fine_index][wall_index] * energy_per_ray
                    
                    # Update wall power variables
                    sub_face.g_a_w[wall_index] = wall_absorbed    # Absorbed power for wall
                    sub_face.e_w[wall_index] = wall_emitted       # Emitted power for wall
                    sub_face.r_w[wall_index] = wall_reflected     # Reflected power for wall
                    sub_face.j_w[wall_index] = wall_emitted + wall_reflected  # Total outgoing power
                    sub_face.g_w[wall_index] = wall_absorbed + wall_reflected  # Total incident power
                    
                    # Wall heat source calculation
                    if sub_face.T_in_w[wall_index] < 0.0
                        # Radiative equilibrium: q_w = 0
                        sub_face.q_w[wall_index] = 0.0
                    else
                        # Boundary element: q_w = absorbed - emitted
                        sub_face.q_w[wall_index] = wall_absorbed - wall_emitted
                    end
                    
                    # Wall temperature calculation
                    if sub_face.T_in_w[wall_index] < 0.0
                        # Radiative equilibrium wall - calculate temperature from emitted power
                        if wall_emitted > 0.0 && sub_face.area[wall_index] > 0.0
                            # For wall: e = epsilon * sigma * T^4 * area
                            # So: T = (e / (epsilon * sigma * area))^(1/4)
                            epsilon_wall = sub_face.epsilon[wall_index]
                            if epsilon_wall > 0.0
                                temp_calc = wall_emitted / (epsilon_wall * STEFAN_BOLTZMANN * sub_face.area[wall_index])
                                if temp_calc > 0.0
                                    sub_face.T_w[wall_index] = temp_calc^(1/4)
                                    if nonzero_temps <= 5  # Print first few
                                        println("Calculated T_w for equilibrium wall ($coarse_index,$fine_index,$wall_index): $(sub_face.T_w[wall_index]) K, emitted = $wall_emitted W")
                                    end
                                end
                            end
                        end
                    else
                        # Prescribed boundary temperature
                        sub_face.T_w[wall_index] = sub_face.T_in_w[wall_index]
                    end
                end
            end
            
            fine_index += 1
        end
    end
    
    println("Temperature calculation summary:")
    println("  Total gas elements: $total_gas_elements")
    println("  Radiative equilibrium elements (T_in_g < 0): $equilibrium_elements")
    println("  Elements with calculated non-zero temperature: $nonzero_temps")
    println("  Maximum temperature: $max_temp K")
    println("  Total gas absorbed power: $total_gas_absorbed_power W")
    println("  Total gas emitted power: $total_gas_emitted_power W")
    println("  Gas power balance: $(total_gas_emitted_power - total_gas_absorbed_power) W")
    println("  gas.kappa = $(gas.kappa)")
end