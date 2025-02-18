using LinearAlgebra

# Vector structure and operations are already defined in Julia
const Point = Vector{Float64} # alias Point to Vector

# Example usage
v1 = [3.0, 4.0, 0.0]
v2 = [1.0, 2.0, 3.0]
println("Vector v1 = ", v1)        # Vector addition
println("Vector v2 = ", v2)        # Vector addition
println("v1 + v2 = ", v1 + v2)        # Vector addition
println("v1 - v2 = ", v1 - v2)        # Vector subtraction
println("Dot product = ", dot(v1, v2))  # Dot product
println(" ⋅  product = ", v1 ⋅ v2)  # Dot product
println("Cross product = ", cross(v1, v2))  # Cross product
println(" × product = ", v1 × v2)  # Cross product
println("Magnitude of v1 = ", norm(v1))  # Magnitude
println("Normalized v1 = ", normalize(v1))  # Normalized vector
