using LinearAlgebra

const Point = Vector # alias Point to Vector
# Define a Ray struct to store the origin and direction
mutable struct Ray
    origin::Point{Float64}
    direction::Vector{Float64}

    function Ray(origin::Vector{Float64}, direction::Vector{Float64})
        new(origin, normalize(direction))  # Normalize direction for consistency
    end
end

# Function to get a point along the ray at parameter t
function at(ray::Ray, t::Float64)
    return ray.origin + t * ray.direction
end

# Example usage of Ray
origin = [0.0, 0.0, 0.0]
direction = [1.0, 0.0, 0.0]
ray = Ray(origin, direction)

# Get the point on the ray at t = 2.0
point = at(ray, 2.0)
if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    println("Ray Origin = ", ray.origin)
    println("Ray Direction = ", ray.direction)
    println("Point on Ray at t = 2.0: ", point)
end