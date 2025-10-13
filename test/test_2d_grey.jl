"""
Test 2D grey participating media ray tracing.
Tests include:
- Square aligned with coordinate system
- Rotated square (45° and other angles)
- Comparison with Crosbie & Schrenker analytical solution
- Mesh refinement independence
"""

println("\n" * "-"^60)
println("Testing 2D Grey Participating Media")
println("-"^60)

# Analytical reference from Crosbie & Schrenker (1984)
# Source function along centerline of 1x1 square with one hot wall
const RELATIVE_TAU_Z = [0.0, 0.00611, 0.02037, 0.04251, 0.07216, 0.10884, 0.15194, 
                        0.20076, 0.25449, 0.31225, 0.37309, 0.43602, 0.50000, 0.56398, 
                        0.62691, 0.68775, 0.74551, 0.79924, 0.84806, 0.89116, 0.92784, 
                        0.95749, 0.97963, 0.99390, 1.00000]

const SOURCE_FUNC_CENTER = [0.6293, 0.6198, 0.6017, 0.5767, 0.5460, 0.5108, 0.4724, 
                           0.4323, 0.3919, 0.3525, 0.3153, 0.2810, 0.2500, 0.2224, 
                           0.1981, 0.1768, 0.1584, 0.1424, 0.1287, 0.1171, 0.1073, 
                           0.0992, 0.0930, 0.0885, 0.0863]

#############################################################################
### HELPER FUNCTIONS #######################################################
#############################################################################

function create_square_domain_2d(; T_hot=1000.0, T_cold=0.0, kappa=1.0, sigma_s=0.0, 
                                  epsilon=1.0, Ndim=11, rotation_angle=0.0)
    """
    Create a 2D square domain for testing.
    rotation_angle in radians rotates the square about its center.
    """
    
    # Base square vertices
    if rotation_angle == 0.0
        vertices = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    else
        # Rotate square around its center (0.5, 0.5)
        cx, cy = 0.5, 0.5
        cos_theta = cos(rotation_angle)
        sin_theta = sin(rotation_angle)
        
        vertices = []
        for (x, y) in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
            # Translate to origin
            x_temp = x - cx
            y_temp = y - cy
            
            # Rotate
            x_rot = x_temp * cos_theta - y_temp * sin_theta
            y_rot = x_temp * sin_theta + y_temp * cos_theta
            
            # Translate back
            push!(vertices, (x_rot + cx, y_rot + cy))
        end
    end
    
    vertices_static = SVector(
        Point2(vertices[1]...),
        Point2(vertices[2]...),
        Point2(vertices[3]...),
        Point2(vertices[4]...)
    )
    
    solidWalls = SVector(true, true, true, true)
    face = PolyFace2D{Float64}(vertices_static, solidWalls, 1, kappa, sigma_s)
    
    # Set boundary conditions
    face.T_in_w = [T_hot, T_cold, T_cold, T_cold]
    face.epsilon = [epsilon, epsilon, epsilon, epsilon]
    face.T_in_g = -1.0  # Unknown gas temperature
    face.q_in_g = 0.0   # Radiative equilibrium
    
    # Create mesh
    mesh = RayTracingMeshOptim([face], [(Ndim, Ndim)])
    
    return mesh
end

function extract_centerline_temperatures(mesh, Ndim)
    """
    Extract temperatures along centerline.
    Returns dimensionless source function along centerline.
    """
    
    # Get all gas temperatures from fine mesh
    all_temps = [fine_face.T_g for fine_face in mesh.fine_mesh[1]]
    
    # Reshape to 2D grid (Ndim x Ndim)
    Tg_matrix = reshape(all_temps, Ndim, Ndim)
    
    # Extract centerline perpendicular to the specified wall
    center_idx = div(Ndim + 1, 2)
    
    # extract centerline
    centerline_temps = Tg_matrix[center_idx, :]
    
    return centerline_temps
end

function dimensionless_source_function(temps, T_hot)
    """Convert temperatures to dimensionless source function"""
    return (temps ./ T_hot).^4
end

"""
Simple linear interpolation function.
Returns a function that interpolates between (x_data, y_data) points.
"""
function line_interpolation(x_data, y_data)
    # Sort data by x values
    sorted_indices = sortperm(x_data)
    x_sorted = x_data[sorted_indices]
    y_sorted = y_data[sorted_indices]
    
    # Return interpolation function
    return function(x)
        if isa(x, Number)
            # Single value
            return _interpolate_single(x, x_sorted, y_sorted)
        else
            # Array of values
            return [_interpolate_single(xi, x_sorted, y_sorted) for xi in x]
        end
    end
end

function _interpolate_single(x, x_data, y_data)
    # Handle extrapolation (constant)
    if x <= x_data[1]
        return y_data[1]
    elseif x >= x_data[end]
        return y_data[end]
    end
    
    # Find bracketing indices
    idx_upper = findfirst(xi -> xi >= x, x_data)
    idx_lower = idx_upper - 1
    
    # Linear interpolation
    x_lower = x_data[idx_lower]
    x_upper = x_data[idx_upper]
    y_lower = y_data[idx_lower]
    y_upper = y_data[idx_upper]
    
    # Interpolate
    t = (x - x_lower) / (x_upper - x_lower)
    return y_lower + t * (y_upper - y_lower)
end

#############################################################################
### TEST 1: ROTATED SQUARE - ALL FOUR WALLS ################################
#############################################################################

@testset "Aligned Square - Single Hot Wall" begin
    T_hot = 1000.0
    T_cold = 0.0
    kappa = 1.0
    sigma_s = 0.0
    Ndim = 11
    epsilon = 1.0
    
    N_rays_total = 1_000_000  # Need more rays to reduce noise like in original
    
    # Create interpolation from analytical solution
    itp_analytical = line_interpolation(RELATIVE_TAU_Z, SOURCE_FUNC_CENTER)
    
    # Sample points (tau values along centerline) - cell centers
    tau_sample = range(1.0/(2*Ndim), 1.0 - 1.0/(2*Ndim), length=Ndim)
    analytical_values = itp_analytical(tau_sample)
    
    # Test each wall as the hot wall
    wall_names = ["Bottom", "Right", "Top", "Left"]
    rotation_angle = 0.0

    for (wall_idx, wall_name) in enumerate(wall_names)
        @testset "$wall_name Wall Hot" begin
            # Create mesh with only this wall hot
            mesh = create_square_domain_2d(T_hot=T_hot, T_cold=T_cold, 
                                          kappa=kappa, sigma_s=sigma_s,
                                          epsilon=epsilon, Ndim=Ndim,
                                          rotation_angle=rotation_angle)
            
            # Set up boundary conditions with only this wall hot
            T_walls = fill(T_cold, 4)
            T_walls[wall_idx] = T_hot
            mesh.coarse_mesh[1].T_in_w = T_walls
            
            # Run ray tracing with exchange factor method
            mesh(N_rays_total; method=:exchange)
            
            # Solve steady state
            steadyStateGrey2D!(mesh, mesh.F_smooth)
            
            # Extract centerline temperatures
            centerline_temps = extract_centerline_temperatures(mesh, Ndim) #, wall_idx)
            
            # Convert to dimensionless source function
            source_func_computed = dimensionless_source_function(centerline_temps, T_hot)
            
            # Compare with analytical solution from Crosbie & Schrenker
            @test isapprox(source_func_computed, analytical_values, rtol=ANALYTICAL_TOLERANCE)
            
            # Energy conservation check
            if !isnothing(mesh.energy_error)
                @test abs(mesh.energy_error) < ENERGY_TOLERANCE # Small relative to total power
            end
        end
        rotation_angle += π/2
    end
end

#############################################################################
### TEST 2: MULTIPLE ROTATION ANGLES #######################################
#############################################################################

@testset "Multiple Rotation Angles" begin
    T_hot = 1000.0
    T_cold = 0.0
    kappa = 1.0
    sigma_s = 0.0
    Ndim = 7  # Smaller for faster testing
    epsilon = 1.0
    N_rays_total = 1_000_000
    
    # Test various rotation angles
    angles = [0.0, π/6, π/4, π/3, π/2, 2π/3]
    
    # Store mean temperatures for comparison
    mean_temps = Float64[]
    
    for angle in angles
        angle_deg = round(rad2deg(angle), digits=1)
        @testset "Rotation $(angle_deg)°" begin
            mesh = create_square_domain_2d(T_hot=T_hot, T_cold=T_cold,
                                          kappa=kappa, sigma_s=sigma_s,
                                          epsilon=epsilon, Ndim=Ndim,
                                          rotation_angle=angle)
            
            # Run ray tracing
            mesh(N_rays_total; method=:exchange)
            
            # Solve steady state
            steadyStateGrey2D!(mesh, mesh.F_smooth)
            
            # Collect statistics
            temps = [fine_face.T_g for fine_face in mesh.fine_mesh[1]]
            T_mean = mean(temps)
            T_min = minimum(temps)
            T_max = maximum(temps)
            
            push!(mean_temps, T_mean)
            
            # Physical bounds
            @test T_min >= 0.0
            @test T_max <= T_hot * (1.0 + ANALYTICAL_TOLERANCE)
            @test T_cold <= T_mean <= T_hot
        end
    end
    
    # Mean temperatures should be relatively consistent across rotations
    @test std(mean_temps) / mean(mean_temps) < ANALYTICAL_TOLERANCE # Less than 5% variation
end

#############################################################################
### TEST 3: MESH REFINEMENT WITH SCALED RAY COUNT ##########################
#############################################################################

@testset "Mesh Refinement with Scaled Rays" begin
    T_hot = 1000.0
    T_cold = 0.0
    kappa = 1.0
    sigma_s = 0.0
    epsilon = 1.0
    
    # For MC methods: need to scale rays with mesh refinement
    # More elements = need more rays to maintain accuracy
    mesh_levels = [5, 7, 9, 11]
    mean_temps = Float64[]
    std_temps = Float64[]
    
    for Ndim in mesh_levels
        @testset "Ndim = $Ndim" begin
            # Scale rays with number of elements to maintain statistical quality
            # Each element needs sufficient rays
            N_elements = Ndim^2 + 4*Ndim
            N_rays_per_element = 2_000
            N_rays_total = N_elements * N_rays_per_element
            
            mesh = create_square_domain_2d(T_hot=T_hot, T_cold=T_cold,
                                          kappa=kappa, sigma_s=sigma_s,
                                          epsilon=epsilon, Ndim=Ndim)
            
            mesh(N_rays_total; method=:exchange)
            steadyStateGrey2D!(mesh, mesh.F_smooth)
            
            temps = [fine_face.T_g for fine_face in mesh.fine_mesh[1]]
            push!(mean_temps, mean(temps))
            push!(std_temps, std(temps))
            
            # Basic sanity checks
            @test all(0.0 .<= temps .<= T_hot * (1.0 + ANALYTICAL_TOLERANCE))
        end
    end
    
    # With properly scaled rays, mean should be relatively constant
    # (not necessarily converging, but consistent within statistical noise)
    @test std(mean_temps) / mean(mean_temps) < ANALYTICAL_TOLERANCE # Mean should be stable
    
end

#############################################################################
### TEST 4: ENERGY CONSERVATION ############################################
#############################################################################

@testset "Energy Conservation" begin
    T_hot = 1000.0
    T_cold = 0.0
    kappa = 1.0
    sigma_s = 0.0
    Ndim = 7
    epsilon = 1.0
    N_rays_total = 1_000_000
    
    mesh = create_square_domain_2d(T_hot=T_hot, T_cold=T_cold,
                                  kappa=kappa, sigma_s=sigma_s,
                                  epsilon=epsilon, Ndim=Ndim)
    
    mesh(N_rays_total; method=:exchange)
    steadyStateGrey2D!(mesh, mesh.F_smooth)
    
    # Check energy error if available
    if !isnothing(mesh.energy_error)
        # Energy error should be small relative to total power
        @test abs(mesh.energy_error) < ENERGY_TOLERANCE # Absolute tolerance in Watts
    end
    
    # Alternative check: sum of all sources and sinks
    total_emission = 0.0
    total_absorption = 0.0
    
    for fine_face in mesh.fine_mesh[1]
        if isa(fine_face.e_g, Number)
            total_emission += fine_face.e_g
            total_absorption += fine_face.g_a_g
        end
    end
    
    # In steady state with radiative equilibrium, these should balance
    rel_error = abs(total_emission - total_absorption) / max(total_emission, 1e-10)
    @test rel_error < ANALYTICAL_TOLERANCE # 5% relative error
end

println("✓ 2D Grey Participating Media tests complete")