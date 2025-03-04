module SourceModule # module begin
using LinearAlgebra
if abspath(PROGRAM_FILE) == @__FILE__
    println("ONLY for DEVELPMENT: Script $(@__FILE__) is running directly!")
    include("./MathModule.jl")  # Include only if necessary
    using .MathModule
    
    # Do something here
else
    # include("./MathModule.jl")  # Include only if necessary
    using Main.MathModule
end

# Define the structure for light sources ===============================================================================



abstract type Source end

struct PointSource <: Source
    # define the point source structure
    frame:: Frame
    wavelength:: Float64 # wavelength of the light [um]
    intensity:: Float64 # intensity of the light source
    distribution:: Function  # Function that returns a sampled direction
end

struct EnsembleSource <: Source
    objects::Vector{Source}
end

mutable struct Ray
    origin::Point       # 3D vector origin position
    direction::Direction # 3D vector direction
    amplitude::Float64  # Amplitude of the ray
    
    # Constructor with default amplitude value
    function Ray(origin::Point, direction::Direction; amplitude::Float64=1.0)
        new(origin, direction, amplitude)
    end
end

function transform_ray_to!(ray::Ray, frame::Frame)
    rotrot = transpose(frame.orientation)
    ray.origin = rotrot*(ray.origin - frame.origin)
    ray.direction = rotrot * ray.direction
end

function transform_ray_from!(ray::Ray, frame::Frame)
    rotrot = frame.orientation
    ray.origin = rotrot*ray.origin + frame.origin
    ray.direction = rotrot * ray.direction
end


function emit_ray(source::Source)::Ray
    direction = source.distribution()
    direction = transform_direction_from(direction, source.frame)
    return Ray(source.frame.origin, direction)
end

function emit_rays(source::Source, num_rays::Int)::Vector{Ray}
    rays = Ray[]
    for _ in 1:num_rays
        push!(rays, emit_ray(source))
    end
    return rays
end

function isotropic_distribution()
    function sample_direction()
        θ = 2π * rand()     # Azimuthal angle (0 to 2π)
        ϕ = acos(2 * rand() - 1)  # Polar angle (cosine-weighted for uniform sphere)
        return [sin(ϕ) * cos(θ), sin(ϕ) * sin(θ), cos(ϕ)]
    end
    return () -> sample_direction()
end

function collimated_distribution(direction::Direction)
    return () -> normalize(direction)  # Always emits in the same direction
end

function cone_distribution(main_direction::Direction, spread_angle::Float64)
    function sample_direction()
        θ = spread_angle * sqrt(rand())  # Cone angle (scaled probability)
        ϕ = 2π * rand()  # Uniform azimuth

        # Local coordinate system
        z = normalize(main_direction)
        x = normalize(cross(z, [0, 1, 0]))  # Ensure perpendicular x
        y = cross(z, x)  # Ensure perpendicular y

        # Convert to Cartesian
        dir = cos(θ) * z + sin(θ) * (cos(ϕ) * x + sin(ϕ) * y)
        return normalize(dir)
    end
    return sample_direction
end

function dipole_distribution(dipole_moment::Direction)
    function sample_direction()
        θ = acos(1 - rand())  # Correct sampling for sin(θ) amplitude distribution
        ϕ = 2π * rand()       # Uniform azimuth

        # Local coordinate system
        z = normalize(dipole_moment)
        x = normalize(cross(z, [0, 1, 0]))  # Ensure perpendicular x
        y = cross(z, x)  # Ensure perpendicular y

        # Convert to Cartesian coordinates
        dir = sin(θ) * (cos(ϕ) * x + sin(ϕ) * y) + cos(θ) * z
        return normalize(dir)
    end
    return sample_direction
end


export Source, PointSource, EnsembleSource, isotropic_distribution, collimated_distribution, cone_distribution, dipole_distribution, 
    Ray, emit_ray, emit_rays, transform_ray_to!, transform_ray_from!

end # module end