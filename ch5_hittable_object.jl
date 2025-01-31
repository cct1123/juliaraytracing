# Import Ray from the previous chapter script
include("ch4.1_ray.jl")  # Assuming it contains the Ray structure and operations

# Define the Sphere struct using Vec3 for the center
mutable struct Sphere
    center::Vec3
    radius::Float64

    function Sphere(center::Vec3, radius::Float64)
        new(center, radius)
    end
end

# Function to compute the intersection of a ray with a sphere using Vec3
function hit(sphere::Sphere, ray::Ray, t_min::Float64, t_max::Float64)
    oc = ray.origin - sphere.center
    a = ray.direction ⋅ ray.direction
    b = oc ⋅ ray.direction
    c = oc ⋅ oc - sphere.radius^2
    discriminant = b^2 - a * c

    if discriminant > 0.0
        temp = (-b - sqrt(discriminant)) / a
        if temp < t_max && temp > t_min
            return true, temp
        end
    end
    return false, 0.0
end

# Example usage of hit function
sphere = Sphere(Vec3(0.0, 0.0, -1.0), 0.5)
ray = Ray(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, -1.0))

hit_result, t = hit(sphere, ray, 0.0, 10.0)
if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    println("Ray hits sphere: ", hit_result)
    println("Intersection at t = ", t)
end
