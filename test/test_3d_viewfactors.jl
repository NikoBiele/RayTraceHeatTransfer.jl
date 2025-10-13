"""
Test 3D view factor calculations against analytical and EES reference values.
Tests include:
- Cube aligned with coordinate system
- Rotated cube (multiple orientations)
- Arbitrary polygon pairs from Narayanaswamy paper
"""

println("\n" * "-"^60)
println("Testing 3D View Factors")
println("-"^60)

#############################################################################
### TEST 1: NARAYANASWAMY EXAMPLES (Arbitrary Polygon Pairs) ###############
#############################################################################

@testset "Narayanaswamy View Factor Examples" begin
    # Reference values from "An analytic expression for radiation view factor 
    # between two arbitrarily oriented planar polygons" by Arvind Narayanaswamy
    
    test_cases = [
        (
            name = "Example 1",
            poly_A = [0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0; 0.0 1.0 0.0],
            poly_B = [0.0 0.0 1.0; 1.0 0.0 1.0; 1.0 1.0 1.0; 0.0 1.0 1.0],
            F_ref = 0.199825
        ),
        (
            name = "Example 2",
            poly_A = [0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0; 0.0 1.0 0.0],
            poly_B = [0.0 0.0 10.0; 1.0 0.0 10.0; 1.0 1.0 10.0; 0.0 1.0 10.0],
            F_ref = 3.16206e-3
        ),
        (
            name = "Example 3",
            poly_A = [0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0; 0.0 1.0 0.0],
            poly_B = [0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 1.0 1.0; 0.0 0.0 1.0],
            F_ref = 0.200044
        ),
        (
            name = "Example 4",
            poly_A = [0.0 0.0 0.0; 0.0 1.0 0.0; 1.0 1.0 0.0],
            poly_B = [1.0 0.0 1.0; 1.0 1.0 1.0; 0.0 1.0 1.0],
            F_ref = 0.099912
        ),
        (
            name = "Example 5a",
            poly_A = [0.0 0.5 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0; 0.0 1.0 0.0],
            poly_B = [2.0 0.5 0.0; 3.0 0.0 0.5; 3.0 2.0 0.5; 2.0 1.5 0.0],
            F_ref = 4.44228e-3
        ),
        (
            name = "Example 5b",
            poly_A = [0.0 0.0 0.0; 0.5 0.0 0.0; 1.0 1.0 0.0; 0.0 1.0 0.0],
            poly_B = [2.0 0.5 0.0; 3.0 0.0 0.5; 3.0 2.0 0.5; 2.0 1.5 0.0],
            F_ref = 3.63699e-3
        ),
        (
            name = "Example 6",
            poly_A = [0.0 0.0 0.0; 1.0 0.0 0.0; 1.0 1.0 0.0],
            poly_B = [2.0 2.0 2.0; 4.0 4.0 4.0; 2.0 3.0 3.0],
            F_ref = 1.06866e-3
        )
    ]
    
    for (i, test) in enumerate(test_cases)
        @testset "$(test.name)" begin
            F_AB, F_BA, area_A, area_B = viewFactor(test.poly_A, test.poly_B)
            
            @test isapprox(F_AB, test.F_ref, atol=VF_TOLERANCE)
            
            # Also test reciprocity: A_a * F_AB = A_b * F_BA
            @test isapprox(area_A * F_AB, area_B * F_BA, rtol=1e-10)
        end
    end
end

#############################################################################
### TEST 2: CUBE VIEW FACTORS (EES Reference) ##############################
#############################################################################

@testset "Cube View Factors - Aligned" begin
    # Create unit cube aligned with coordinate system
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
        1 2 4 3;  # Front face
        5 6 8 7;  # Back face
        1 5 7 3;  # Bottom face
        2 6 8 4;  # Top face
        3 4 8 7;  # Left face
        1 2 6 5   # Right face
    ]
    
    # Accurately converged view factors from EES (for 1x1x1 cube)
    F_EES = [
        0.000000000000000000  0.199824895698387383  0.200043776075403154  0.200043776075403154  0.200043776075403154  0.200043776075403154;
        0.199824895698387383  0.000000000000000000  0.200043776075403154  0.200043776075403154  0.200043776075403154  0.200043776075403154;
        0.200043776075403154  0.200043776075403154  0.000000000000000000  0.199824895698387383  0.200043776075403154  0.200043776075403154;
        0.200043776075403154  0.200043776075403154  0.199824895698387383  0.000000000000000000  0.200043776075403154  0.200043776075403154;
        0.200043776075403154  0.200043776075403154  0.200043776075403154  0.200043776075403154  0.000000000000000000  0.199824895698387383;
        0.200043776075403154  0.200043776075403154  0.200043776075403154  0.200043776075403154  0.199824895698387383  0.000000000000000000
    ]
    
    Ndim = 1
    epsilon = ones(size(faces, 1))
    q_in_w = zeros(size(faces, 1))
    T_in_w = -ones(size(faces, 1))  # All unknown
    
    # Create domain and compute view factors
    domain3D = Domain3D_faces(points, faces, Ndim, q_in_w, T_in_w, epsilon)
    
    # Test against EES reference
    @test maximum(abs.(domain3D.F - F_EES)) < VF_TOLERANCE
    
    # Test reciprocity for all face pairs
    for i in 1:6, j in 1:6
        if i != j
            area_i = sum([sf.area for sf in domain3D.facesMesh[i].subFaces])
            area_j = sum([sf.area for sf in domain3D.facesMesh[j].subFaces])
            reciprocity_error = abs(area_i * domain3D.F[i,j] - area_j * domain3D.F[j,i])
            @test reciprocity_error < 1e-10
        end
    end
    
    # Test summation rule: sum of view factors from each face should be ~1
    for i in 1:6
        @test isapprox(sum(domain3D.F[i,:]), 1.0, atol=1e-10)
    end
end

#############################################################################
### TEST 3: ROTATED CUBE VIEW FACTORS ######################################
#############################################################################

function rotate_points(points, axis, angle)
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

@testset "Cube View Factors - Rotations" begin
    # Base cube
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
    
    # Reference view factors from aligned cube
    F_EES = [
        0.000000000000000000  0.199824895698387383  0.200043776075403154  0.200043776075403154  0.200043776075403154  0.200043776075403154;
        0.199824895698387383  0.000000000000000000  0.200043776075403154  0.200043776075403154  0.200043776075403154  0.200043776075403154;
        0.200043776075403154  0.200043776075403154  0.000000000000000000  0.199824895698387383  0.200043776075403154  0.200043776075403154;
        0.200043776075403154  0.200043776075403154  0.199824895698387383  0.000000000000000000  0.200043776075403154  0.200043776075403154;
        0.200043776075403154  0.200043776075403154  0.200043776075403154  0.200043776075403154  0.000000000000000000  0.199824895698387383;
        0.200043776075403154  0.200043776075403154  0.200043776075403154  0.200043776075403154  0.199824895698387383  0.000000000000000000
    ]
    
    # Test various rotations
    rotations = [
        (:x, π/4, "45° around X"),
        (:y, π/4, "45° around Y"),
        (:z, π/4, "45° around Z"),
        (:x, π/3, "60° around X"),
        (:y, π/6, "30° around Y"),
        (:z, π/2, "90° around Z")
    ]
    
    Ndim = 1
    epsilon = ones(6)
    q_in_w = zeros(6)
    T_in_w = -ones(6)
    
    for (axis, angle, desc) in rotations
        @testset "Rotation: $desc" begin
            # Rotate the cube
            points_rotated = rotate_points(points_base, axis, angle)
            
            # Create domain with rotated geometry
            domain3D = Domain3D_faces(points_rotated, faces, Ndim, q_in_w, T_in_w, epsilon)
            
            # Extract unique view factor values (excluding self-view)
            F_unique = Float64[]
            for i in 1:6, j in i+1:6
                push!(F_unique, domain3D.F[i,j])
            end
            sort!(F_unique)
            
            F_EES_unique = Float64[]
            for i in 1:6, j in i+1:6
                push!(F_EES_unique, F_EES[i,j])
            end
            sort!(F_EES_unique)
            
            # The sorted unique values should match
            @test isapprox(F_unique, F_EES_unique, atol=VF_TOLERANCE)
            
            # Test reciprocity still holds
            for i in 1:6, j in 1:6
                if i != j
                    area_i = sum([sf.area for sf in domain3D.facesMesh[i].subFaces])
                    area_j = sum([sf.area for sf in domain3D.facesMesh[j].subFaces])
                    @test isapprox(area_i * domain3D.F[i,j], area_j * domain3D.F[j,i], rtol=1e-10)
                end
            end
            
            # Test summation rule
            for i in 1:6
                @test isapprox(sum(domain3D.F[i,:]), 1.0, atol=1e-10)
            end
        end
    end
end

println("✓ 3D View Factor tests complete")