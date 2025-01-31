using LinearAlgebra

# Import Vec3 and Ray from ch3.1_ray_implementation.jl
include("ch4.1j_ray.jl")

# Sphere structure to define a sphere and its radius
struct Sphere
    center::Vector
    radius::Float64

    function Sphere(center::Vector, radius::Float64)
        new(center, radius)
    end
end

# Function to compute intersection with a sphere
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

# Example usage of Sphere and hit function
sphere = Sphere([0.0, 0.0, -1.0], 0.5)
ray = Ray([0.0, 0.0, 0.0], [0.0, 0.0, -1.0])

# Check if the ray hits the sphere
hit_result, t = hit(sphere, ray, 0.0, 10.0)

if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    println("Ray hits sphere: ", hit_result)
    println("Intersection at t = ", t)
end
