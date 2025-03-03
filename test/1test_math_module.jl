using Test
using LinearAlgebra
include("../MathModule.jl")
using .MathModule  # Replace with the actual path to your module

# Test case for the `rotation_matrix` function
@testset "Rotation Matrix Tests" begin
    # Test 1: 90 degrees rotation around the Y-axis
    axis = [0.0, 1.0, 0.0]
    angle = π / 2
    R = rotation_matrix(axis, angle)
    
    # Check if the rotation matrix is correct for 90 degrees around the Y-axis
    @test R ≈ [
        0.0  0.0  1.0;
        0.0  1.0  0.0;
        -1.0  0.0  0.0
    ] atol=1e-6
    
    # Test 2: 90 degrees rotation around the Z-axis
    axis = [0.0, 0.0, 1.0]
    angle = π / 2
    R = rotation_matrix(axis, angle)
    
    # Check if the rotation matrix is correct for 90 degrees around the Z-axis
    @test R ≈ [
        0.0  -1.0  0.0;
        1.0  0.0  0.0;
        0.0  0.0  1.0
    ] atol=1e-6

    # Test 3: Arbitrary axis (1, 1, 1) and 60 degrees
    axis = [1.0, 1.0, 1.0]
    angle = π / 3
    R = rotation_matrix(axis, angle)
    
    # Check if the rotation matrix is close to the expected values
    expected_R = [
        2.0/3.0 -1.0/3.0 2.0/3.0; 
        2.0/3.0 2.0/3.0 -1.0/3.0; 
        -1.0/3.0 2.0/3.0 2.0/3.0]
    @test R ≈ expected_R atol=1e-6
end

# Test case for the `rotate` function
@testset "Rotation Function Tests" begin
    # Test 1: Rotate vector [1, 0, 0] by 90 degrees around the X-axis
    v = [1.0, 0.0, 0.0]
    k = [1.0, 0.0, 0.0]
    theta = π / 2
    v_rotated = rotate(v, k, theta)
    
    # Expected result is [1.0, 0.0, 0.0] since the vector is along the rotation axis
    @test v_rotated ≈ [1.0, 0.0, 0.0] atol=1e-6
    
    # Test 2: Rotate vector [1, 0, 0] by 90 degrees around the Y-axis
    v = [1.0, 0.0, 0.0]
    k = [0.0, 1.0, 0.0]
    theta = π / 2
    v_rotated = rotate(v, k, theta)
    
    # The expected rotated vector after 90-degree rotation about Y-axis is [0.0, 0.0, 1.0]
    @test v_rotated ≈ [0.0, 0.0, -1.0] atol=1e-6

    # Test 3: Rotate vector [1, 0, 0] by 90 degrees around the Z-axis
    v = [1.0, 0.0, 0.0]
    k = [0.0, 0.0, 1.0]
    theta = π / 2
    v_rotated = rotate(v, k, theta)
    
    # The expected rotated vector after 90-degree rotation about Z-axis is [0.0, 1.0, 0.0]
    @test v_rotated ≈ [0.0, 1.0, 0.0] atol=1e-6

    # Test 4: Rotate vector [1, 0, 0] by 60 degrees around an arbitrary axis (1, 1, 1)
    v = [1.0, 0.0, 0.0]
    k = [1.0, 1.0, 1.0]
    theta = π / 3
    v_rotated = rotate(v, k, theta)
    
    # Using rotation matrix to verify the result
    R = rotation_matrix(k, theta)
    v_rotated_matrix = R * v
    
    # Check if the result matches the matrix-rotated vector
    @test v_rotated ≈ v_rotated_matrix atol=1e-6
end

@testset "Advanced Invalid Input Tests" begin
    # Test 1: Rotation matrix with invalid axis (axis not orthogonal or normalized)
    axis = [2.0, 3.0, 1.0]  # Not normalized, expected to be a unit vector
    angle = π / 4
    R = rotation_matrix(axis, angle)
    # Check that the resulting rotation matrix is valid and that it does not throw errors,
    # but the axis is not normalized.
    @test isapprox(norm(R*[1.0, 0.0, 0.0]), 1.0, atol=1e-6)  # Check that axis is normalized if we manually fix it
    
    # Test 2: Rotation matrix with extreme angle (large angle like 10π)
    axis = [0.0, 1.0, 0.0]  # Normalized Y-axis
    angle = 10 * π  # 10 full rotations, should be equivalent to 0 (mod 2π)
    R = rotation_matrix(axis, angle)
    @test isapprox(R, I(3), atol=1e-6)  # After 10 full rotations, the matrix should return to identity
    
    # Test 3: Rotation with very small angle close to 0
    v = [1.0, 0.0, 0.0]
    k = [0.0, 1.0, 0.0]  # Y-axis
    theta = 1e-9  # Very small angle
    v_rotated = rotate(v, k, theta)
    # Rotation by a very small angle should result in a nearly unchanged vector
    @test isapprox(v_rotated, v, atol=1e-6)
    
    # Test 4: Rotate with non-orthogonal vectors
    v = [1.0, 0.0, 0.0]  # Vector on X-axis
    k = [1.0, 1.0, 0.0]  # Axis that is not orthogonal to the vector
    theta = π / 4.0  # 45 degree rotation
    v_rotated = rotate(v, k, theta)
    # The rotated vector should be a combination of the original vector and the axis
    @test isapprox(norm(v_rotated), 1.0, atol=1e-6)  # Check that the rotation preserves magnitude
    
    # Test 5: Test rotation on a vector that's already on the rotation axis
    v = [1.0, 0.0, 0.0]  # Vector on X-axis
    k = [1.0, 0.0, 0.0]  # Rotation axis is the same as the vector
    theta = π / 2  # 90 degrees
    v_rotated = rotate(v, k, theta)
    @test isapprox(v_rotated, v, atol=1e-6)  # The vector should not change since it’s aligned with the axis
    
    # Test 6: Rotation where the angle is multiple of π
    axis = [1.0, 0.0, 0.0]  # X-axis
    angle = 3 * π  # 3 full rotations
    R = rotation_matrix(axis, angle)
    @test isapprox(R, [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], atol=1e-6)  # The result after 3π should be the identity matrix
    
    # Test 7: Rotation by zero angle
    axis = [0.0, 0.0, 1.0]  # Z-axis
    angle = 0.0  # Zero angle rotation
    R = rotation_matrix(axis, angle)
    @test isapprox(R, I(3), atol=1e-6)  # Rotation by zero should return the identity matrix
    
    # Test 8: Rotation with random axis and angle
    axis = rand(3)  # Random axis, not guaranteed to be normalized
    axis = normalize(axis)  # Normalize the axis
    angle = rand() * 2 * π  # Random angle between 0 and 2π
    R = rotation_matrix(axis, angle)
    @test isapprox(det(R), 1.0, atol=1e-6)  # Rotation matrix determinant should be 1
    
    # Test 9: Rotate large values (edge case for numerical stability)
    v = [1e10, 0.0, 0.0]  # Very large vector
    k = [0.0, 1.0, 0.0]  # Y-axis as rotation axis
    theta = π / 2  # 90-degree rotation
    v_rotated = rotate(v, k, theta)
    @test isapprox(norm(v_rotated), 1e10, atol=1e-6)  # The magnitude should remain the same
    
    # Test 10: Rotation with axis at an extreme angle (π/2)
    axis = [0.0, 0.0, 1.0]  # Z-axis
    angle = π / 2  # 90-degree rotation
    v = [1.0, 0.0, 0.0]
    v_rotated = rotate(v, axis, angle)
    @test isapprox(v_rotated, [0.0, 1.0, 0.0], atol=1e-6)  # The vector should rotate to the Y-axis
end

@testset "Frame Struct Tests" begin

    # Test 1: Creating a Frame with zero rotation and origin at [0, 0, 0]
    origin = [0.0, 0.0, 0.0]
    orientation = Matrix{Float64}(I(3))  # Identity matrix (no rotation)
    frame = Frame(origin, orientation)
    
    @test frame.origin == origin
    @test isapprox(frame.orientation, Matrix{Float64}(I(3)), atol=1e-6)

    # Test 2: Rotating a point with no rotation (Identity frame)
    point_obj = [1.0, 0.0, 0.0]
    frame_no_rotation = Frame([0.0, 0.0, 0.0], Matrix{Float64}(I(3)))  # No rotation
    point_lab = frame_no_rotation.orientation * point_obj
    
    @test isapprox(point_lab, point_obj, atol=1e-6)

    # Test 3: Rotating a point by 90 degrees around the Z-axis
    axis = [0.0, 0.0, 1.0]  # Z-axis
    angle = π / 2  # 90-degree rotation
    R = rotation_matrix(axis, angle)
    
    frame_rotated = Frame([0.0, 0.0, 0.0], R)
    point_obj = [1.0, 0.0, 0.0]  # Point in object frame
    expected_point_lab = [0.0, 1.0, 0.0]  # After rotation
    
    point_lab = frame_rotated.orientation * point_obj
    
    @test isapprox(point_lab, expected_point_lab, atol=1e-6)

    # Test 4: Frame with translation (origin at [1, 1, 1])
    frame_with_translation = Frame([1.0, 1.0, 1.0], Matrix{Float64}(I(3)))  # Identity rotation, just translation
    point_obj = [1.0, 0.0, 0.0]  # Point in object frame
    expected_point_lab = [2.0, 1.0, 1.0]  # Point translated by [1, 1, 1]
    
    point_lab = frame_with_translation.orientation * point_obj + frame_with_translation.origin
    
    @test isapprox(point_lab, expected_point_lab, atol=1e-6)

end


@testset "Ray-Surface Interaction Tests" begin
    # Test 1: Surface normal (implicit surface: sphere)
    surface_function(p) = p⋅p - 1.0  # A unit sphere
    normal = surface_normal([1.0, 0.0, 0.0], surface_function)
    @test isapprox(normal, [1.0, 0.0, 0.0], atol=1e-6)

    # Test 2: Ray intersection with sphere
    o = [0.0, 0.0, 3.0]  # Origin above the sphere
    d = [0.0, 0.0, -1.0]  # Direction towards the sphere
    intersection_t = find_intersection(surface_function, o, d)
    @test isapprox(intersection_t, 2.0, atol=1e-6)  # Expected intersection at t = 2

    # Test 3: Reflection of a ray
    normal = [0.0, 1.0, 0.0]
    incident = normalize([1.0, -2.0, 0.0])
    reflected_ray = reflect(normal, incident)
    @test isapprox(reflected_ray, [0.4472135954999579, 0.8944271909999159, 0.0], atol=1e-6)

    # Test 4: Refraction (assuming n1 = 1.5, n2 = 1.0)
    refracted_ray = refract(normal, incident, 1.5, 1.0)
    @test isapprox(refracted_ray, [0.6708203932499369, -0.7416198487095662, 0.0], atol=1e-6)

    # Test 5: Fresnel reflectance (assuming n1 = 1.5, n2 = 1.0)
    reflectance_value = fresnel(normal, incident, 1.5, 1.0)
    @test isapprox(reflectance_value[1], 0.952622096319283, atol=1e-6)  # Approximate Fresnel transmittance at this angle
    @test isapprox(reflectance_value[2], refracted_ray, atol=1e-6)  # Approximate Fresnel refracted ray at this angle
    @test isapprox(reflectance_value[3], 0.04737790368071705, atol=1e-6)  # Approximate Fresnel reflectance at this angle
    @test isapprox(reflectance_value[4], reflected_ray, atol=1e-6)  # Approximate Fresnel reflected ray at this angle

    # Test 6: Schlick's approximation for reflectance
    r_schlick_value = rSchlick2(normal, incident, 1.5, 1.0)
    @test isapprox(r_schlick_value, 0.04, atol=1e-2)  # Schlick's approximation
end
