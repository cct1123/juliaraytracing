# Import Vec3 from the previous chapter script
include("ch3_vec3.jl")  # Assuming it contains the Vec3 structure and operations

# Ray struct to store origin and direction using Vec3
mutable struct Ray
    origin::Point3
    direction::Vec3

    function Ray(origin::Vec3, direction::Vec3)
        new(origin, normalize(direction))  # Normalize direction for consistency
    end
end

# Function to get a point along the ray at parameter t
function at(ray::Ray, t::Float64)
    return ray.origin + t * ray.direction
end

# Example usage of Ray
origin = Vec3(0.0, 0.0, 0.0)
direction = Vec3(1.0, 0.0, 0.0)
ray = Ray(origin, direction)


# Get the point on the ray at t = 2.0
point = at(ray, 2.0)


if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    println("Ray Origin = ", ray.origin)
    println("Ray Direction = ", ray.direction)

    println("Point on Ray at t = 2.0: ", point)
end