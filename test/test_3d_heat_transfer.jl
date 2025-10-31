"""
Test 3D heat transfer calculations with transparent surfaces.
Tests include:
- Isothermal enclosure (all walls at same temperature)
- Energy conservation
- Rotational invariance of solutions
"""

println("\n" * "-"^60)
println("Testing 3D Heat Transfer")
println("-"^60)

#############################################################################
### TEST 1: ISOTHERMAL CUBE (All walls at same temperature) ###############
#############################################################################

@testset "Isothermal Cube" begin
    points = [
        0.0 0.0 0.0;
        0.0 0.0 1.0;
        0.0 1.0 0.0; 
        0.0 1.0 1.0; 
        1.0 0.0 0.0; 
        1.0 0.0 1.0; 
        1.0 1.0 0.0; 
        1.0 1.0 1.0
    ]
    
    faces = [
        1 2 4 3;
        5 6 8 7;
        1 5 7 3;
        2 6 8 4;
        3 4 8 7;
        1 2 6 5
    ]
    
    Ndim = 5
    T_iso = 1000.0  # K
    epsilon = ones(6)
    
    # All walls at same temperature, zero heat flux
    q_in_w = zeros(6)
    T_in_w = fill(T_iso, 6)
    
    domain3D = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w, epsilon)
    steadyStateGrey3D!(domain3D, domain3D.F)
    
    # Extract temperatures from all subfaces
    for (i, superface) in enumerate(domain3D.facesMesh)
        for subface in superface.subFaces
            # All surfaces should remain at isothermal temperature
            @test isapprox(subface.T_w, T_iso, atol=TEMP_TOLERANCE)
            
            # Net heat flux should be zero (within numerical precision)
            @test abs(subface.q_w) < ENERGY_TOLERANCE
        end
    end
end

#############################################################################
### TEST 2: TWO HOT WALLS (Opposite walls at different temperatures) ######
#############################################################################

@testset "Two Hot Walls - Opposite" begin
    points = [
        0.0 0.0 0.0;
        0.0 0.0 1.0;
        0.0 1.0 0.0; 
        0.0 1.0 1.0; 
        1.0 0.0 0.0; 
        1.0 0.0 1.0; 
        1.0 1.0 0.0; 
        1.0 1.0 1.0
    ]
    
    faces = [
        1 2 4 3;  # Front (face 1)
        5 6 8 7;  # Back (face 2)
        1 5 7 3;  # Bottom
        2 6 8 4;  # Top
        3 4 8 7;  # Left
        1 2 6 5   # Right
    ]
    
    Ndim = 5
    T_hot = 1000.0  # K
    T_cold = 500.0  # K
    epsilon = ones(6)
    
    # Front and back walls at different temperatures
    T_in_w = [T_hot, T_cold, -1.0, -1.0, -1.0, -1.0]
    q_in_w = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    domain3D = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w, epsilon)
    steadyStateGrey3D!(domain3D, domain3D.F)
    
    # Check that specified temperatures are maintained
    avg_T_face1 = mean([sf.T_w for sf in domain3D.facesMesh[1].subFaces])
    avg_T_face2 = mean([sf.T_w for sf in domain3D.facesMesh[2].subFaces])
    
    @test isapprox(avg_T_face1, T_hot, atol=TEMP_TOLERANCE)
    @test isapprox(avg_T_face2, T_cold, atol=TEMP_TOLERANCE)
    
    # Other walls should have temperatures between hot and cold
    for i in 3:6
        avg_T = mean([sf.T_w for sf in domain3D.facesMesh[i].subFaces])
        @test T_cold < avg_T < T_hot
    end
    
    # Energy conservation: total heat in should equal total heat out
    total_q = sum(sf.q_w for i in 1:6 for sf in domain3D.facesMesh[i].subFaces)
    @test abs(total_q) < ENERGY_TOLERANCE  # Should sum to nearly zero
end

#############################################################################
### TEST 3: ENERGY CONSERVATION WITH SPECIFIED HEAT FLUX ###################
#############################################################################

@testset "Energy Conservation - Heat Flux BCs" begin
    points = [
        0.0 0.0 0.0;
        0.0 0.0 1.0;
        0.0 1.0 0.0; 
        0.0 1.0 1.0; 
        1.0 0.0 0.0; 
        1.0 0.0 1.0; 
        1.0 1.0 0.0; 
        1.0 1.0 1.0
    ]
    
    faces = [
        1 2 4 3;
        5 6 8 7;
        1 5 7 3;
        2 6 8 4;
        3 4 8 7;
        1 2 6 5
    ]
    
    Ndim = 5
    epsilon = ones(6)
    
    # Specify heat flux on two walls, solve for others
    q_heating = 1000.0  # W
    T_in_w = [1000.0, 500.0, -1.0, -1.0, -1.0, -1.0]
    q_in_w = [0.0, 0.0, q_heating, q_heating, 0.0, 0.0]
    
    domain3D = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w, epsilon)
    steadyStateGrey3D!(domain3D, domain3D.F)
    
    # Calculate total energy balance
    total_q_in = 0.0
    total_q_out = 0.0
    
    for (i, superface) in enumerate(domain3D.facesMesh)
        face_q_total = sum(sf.q_w for sf in superface.subFaces)
        
        if face_q_total > 0
            total_q_in += face_q_total
        else
            total_q_out += abs(face_q_total)
        end
    end
    
    # Energy in should approximately equal energy out
    rel_error = abs(total_q_in - total_q_out) / max(total_q_in, total_q_out)
    @test rel_error < ENERGY_TOLERANCE
end

#############################################################################
### TEST 4: ROTATIONAL INVARIANCE ##########################################
#############################################################################

function rotate_points_3d(points, axis, angle)
    """Rotate points around given axis by angle (radians)"""
    if axis == :x
        R = [1.0 0.0 0.0;
             0.0 cos(angle) -sin(angle);
             0.0 sin(angle) cos(angle)]
    elseif axis == :y
        R = [cos(angle) 0.0 sin(angle);
             0.0 1.0 0.0;
             -sin(angle) 0.0 cos(angle)]
    elseif axis == :z
        R = [cos(angle) -sin(angle) 0.0;
             sin(angle) cos(angle) 0.0;
             0.0 0.0 1.0]
    else
        error("Unknown axis: $axis")
    end
    
    return points * R'
end

@testset "Rotational Invariance of Solutions" begin
    # Base configuration
    points_base = [
        0.0 0.0 0.0;
        0.0 0.0 1.0;
        0.0 1.0 0.0; 
        0.0 1.0 1.0; 
        1.0 0.0 0.0; 
        1.0 0.0 1.0; 
        1.0 1.0 0.0; 
        1.0 1.0 1.0
    ]
    
    faces = [
        1 2 4 3;
        5 6 8 7;
        1 5 7 3;
        2 6 8 4;
        3 4 8 7;
        1 2 6 5
    ]
    
    Ndim = 4  # Smaller mesh for faster testing
    epsilon = ones(6)
    T_hot = 1000.0
    T_cold = 500.0
    
    # Solve base case: front wall hot, back wall cold
    T_in_w = [T_hot, T_cold, -1.0, -1.0, -1.0, -1.0]
    q_in_w = zeros(6)
    
    domain_base = ViewFactorDomain3D(points_base, faces, Ndim, q_in_w, T_in_w, epsilon)
    steadyStateGrey3D!(domain_base, domain_base.F)
    
    # Get temperature statistics from base case
    base_T_min = minimum(sf.T_w for i in 1:6 for sf in domain_base.facesMesh[i].subFaces)
    base_T_max = maximum(sf.T_w for i in 1:6 for sf in domain_base.facesMesh[i].subFaces)
    base_T_mean = mean(sf.T_w for i in 1:6 for sf in domain_base.facesMesh[i].subFaces)
    
    # Test rotations - solution statistics should be similar
    rotations = [
        (:z, π/4, "45° around Z"),
        (:x, π/6, "30° around X"),
        (:y, π/3, "60° around Y")
    ]
    
    for (axis, angle, desc) in rotations
        @testset "Rotation: $desc" begin
            points_rot = rotate_points_3d(points_base, axis, angle)
            
            domain_rot = ViewFactorDomain3D(points_rot, faces, Ndim, q_in_w, T_in_w, epsilon)
            steadyStateGrey3D!(domain_rot, domain_rot.F)
            
            # Temperature statistics should be invariant
            rot_T_min = minimum(sf.T_w for i in 1:6 for sf in domain_rot.facesMesh[i].subFaces)
            rot_T_max = maximum(sf.T_w for i in 1:6 for sf in domain_rot.facesMesh[i].subFaces)
            rot_T_mean = mean(sf.T_w for i in 1:6 for sf in domain_rot.facesMesh[i].subFaces)
            
            @test isapprox(rot_T_min, base_T_min, rtol=0.01)
            @test isapprox(rot_T_max, base_T_max, rtol=0.01)
            @test isapprox(rot_T_mean, base_T_mean, rtol=0.01)
            
            # Energy conservation should still hold
            total_q = sum(sf.q_w for i in 1:6 for sf in domain_rot.facesMesh[i].subFaces)
            @test abs(total_q) < ENERGY_TOLERANCE
        end
    end
end

#############################################################################
### TEST 5: GREY SURFACE PROPERTIES #######################################
#############################################################################

@testset "Grey Surface Properties" begin
    points = [
        0.0 0.0 0.0;
        0.0 0.0 1.0;
        0.0 1.0 0.0; 
        0.0 1.0 1.0; 
        1.0 0.0 0.0; 
        1.0 0.0 1.0; 
        1.0 1.0 0.0; 
        1.0 1.0 1.0
    ]
    
    faces = [
        1 2 4 3;
        5 6 8 7;
        1 5 7 3;
        2 6 8 4;
        3 4 8 7;
        1 2 6 5
    ]
    
    Ndim = 4
    
    # Test with non-black surfaces
    epsilon = [0.8, 0.8, 0.6, 0.6, 0.9, 0.9]
    T_in_w = [1000.0, 500.0, -1.0, -1.0, -1.0, -1.0]
    q_in_w = zeros(6)
    
    domain3D = ViewFactorDomain3D(points, faces, Ndim, q_in_w, T_in_w, epsilon)
    steadyStateGrey3D!(domain3D, domain3D.F)
    
    # Solution should still exist and be physically reasonable
    for (i, superface) in enumerate(domain3D.facesMesh)
        for subface in superface.subFaces
            # Temperature should be between the specified extremes
            @test 400.0 < subface.T_w < 1100.0
            
            # Heat flux should be finite
            @test isfinite(subface.q_w)
        end
    end
    
    # Energy conservation
    total_q = sum(sf.q_w for i in 1:6 for sf in domain3D.facesMesh[i].subFaces)
    @test abs(total_q) < ENERGY_TOLERANCE
end

println("✓ 3D Heat Transfer tests complete")