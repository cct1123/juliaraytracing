function rotation_matrix(axis::Vector{Float64}, angle::Float64)::Matrix{Float64}
    # Ensure the axis is a unit vector
    axis = normalize(axis)
    
    # Extract components of the axis
    kx, ky, kz = axis
    
    # Compute trigonometric values
    cos_theta = cos(angle)
    sin_theta = sin(angle)
    one_minus_cos = 1 - cos_theta
    
    # Construct the rotation matrix using Rodrigues' formula
    R = [
        cos_theta + kx^2 * one_minus_cos         kx * ky * one_minus_cos - kz * sin_theta   kx * kz * one_minus_cos + ky * sin_theta;
        ky * kx * one_minus_cos + kz * sin_theta  cos_theta + ky^2 * one_minus_cos          ky * kz * one_minus_cos - kx * sin_theta;
        kz * kx * one_minus_cos - ky * sin_theta  kz * ky * one_minus_cos + kx * sin_theta   cos_theta + kz^2 * one_minus_cos
    ]
    
    return R
end

function rotate(v::Vector{Float64}, k::Vector{Float64}, theta::Float64)::Vector{Float64}
    # Normalize the rotation axis vector k
    k = normalize(k)
    
    # Cross product of k and v
    k_cross_v = cross(k, v)
    
    # Dot product of k and v
    k_dot_v = dot(k, v)
    
    # Apply Rodrigues' rotation formula
    v_rotated = v * cos(theta) + k_cross_v * sin(theta) + k * k_dot_v * (1 - cos(theta))
    
    return v_rotated
end

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    # Example Usage for rotation_matrix------------------------------------
    axis = [0.0, 1.0, 0.0]  # Y-axis
    angle = π / 2           # 90 degrees in radians

    R = rotation_matrix(axis, angle)
    println("Rotation Matrix (Y-axis, 90°):")
    println(R)

    axis = [0.0, 0.0, 1.0]  # Z-axis
    angle = π / 2           # 90 degrees in radians

    R = rotation_matrix(axis, angle)
    println("Rotation Matrix (Z-axis, 90°):")
    println(R)

    axis = [1.0, 1.0, 1.0]  # Arbitrary axis
    angle = π / 3           # 60 degrees in radians

    R = rotation_matrix(axis, angle)
    println("Rotation Matrix (Arbitrary axis, 60°):")
    println(R)

    # Verify orthogonality
    println("\nOrthogonality Check (R^T * R):")
    println(R' * R)

    # Verify determinant
    println("\nDeterminant of R:")
    println(det(R))

    # Example Usage for rotate-----------------------------------------
    v = [1.0, 0.0, 0.0]  # Vector to rotate
    k = [1.0, 0.0, 0.0]  # Rotation axis (X-axis)
    theta = π / 2        # 90 degrees in radians
    
    v_rotated = rotate(v, k, theta)
    println("Rotated Vector (X-axis, 90°):")
    println(v_rotated)
    
    v = [1.0, 0.0, 0.0]  # Vector to rotate
    k = [0.0, 1.0, 0.0]  # Rotation axis (Y-axis)
    theta = π / 2        # 90 degrees in radians
    
    v_rotated = rotate(v, k, theta)
    println("Rotated Vector (Y-axis, 90°):")
    println(v_rotated)
    
    v = [1.0, 0.0, 0.0]  # Vector to rotate
    k = [0.0, 0.0, 1.0]  # Rotation axis (Z-axis)
    theta = π / 2        # 90 degrees in radians
    
    v_rotated = rotate(v, k, theta)
    println("Rotated Vector (Z-axis, 90°):")
    println(v_rotated)
    
    v = [1.0, 0.0, 0.0]  # Vector to rotate
    k = [1.0, 1.0, 1.0]  # Rotation axis (arbitrary)
    theta = π / 3        # 60 degrees in radians
    
    # Rotate using the function
    v_rotated = rotate(v, k, theta)
    println("Rotated Vector (Arbitrary axis, 60°):")
    println(v_rotated)
    
    # Verify using the rotation matrix
    R = rotation_matrix(k, theta)
    v_rotated_matrix = R * v
    println("\nRotated Vector (Using Rotation Matrix):")
    println(v_rotated_matrix)


    
end
