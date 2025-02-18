# rtweekend.jl
module RTWeekend
    # Include necessary files
    include("vec3.jl")  # Include Vec3 and related functions
    include("ray.jl")   # Include Ray struct and related functions
    include("color.jl") # Include color utilities
    include("interval.jl")
    include("ppm_to_jpg.jl")
    # Constants
    const infinity = Inf64  # Represents an infinitely large value
    const pi = π            # Julia already has π built-in, but we can alias it for clarity

    # Utility Functions
    degrees_to_radians(degrees) = degrees * pi / 180.0
    function random_double()
        rand()
    end
    
    function random_double(min::Float64, max::Float64)
        min + (max - min) * random_double()
    end
    
    function random_vec3(min::Float64, max::Float64)
        return Vec3(random_double(min, max), random_double(min, max), random_double(min, max))
    end
    
    function length_squared(v::Vec3)
        return norm(v)^2
    end
    
    function unit_vector(v::Vec3)
        return normalize(v)
    end
    
    function near_zero(v::Vec3)
        s = 1e-8
        return abs(v.x) < s && abs(v.y) < s && abs(v.z) < s
    end
    
    function random_unit_vector()
        while true
            p = random_vec3(-1.0, 1.0)
            len_sq = length_squared(p)
            if 1e-160 < len_sq <= 1.0
                return unit_vector(p)
            end
        end
    end
    
    function random_on_hemisphere(normal::Vec3)
        on_unit_sphere = random_unit_vector()
        if dot(on_unit_sphere, normal) > 0.0
            return on_unit_sphere
        else
            return -on_unit_sphere
        end
    end

    function random_in_unit_disk()
        while true
            p = Vec3(random_double(-1.0, 1.0), random_double(-1.0, 1.0), 0.0)
            if length_squared(p) < 1.0
                return p
            end
        end
    end

    
    function my_clamp(x, min_val, max_val)
        if x < min_val
            return min_val
        elseif x > max_val
            return max_val
        else
            return x
        end
    end
    
    # Automatically export all symbols defined in this module
    for sym in names(@__MODULE__; all=true)
        # Skip unwanted symbols
        if !(sym in [
            :RTWeekend,               # Module name
            Symbol("@__MODULE__"),    # Internal symbol
            :eval,                    # Built-in function
            :include,                 # Built-in function
            :names,                   # Built-in function
            :rand,                    # Built-in function
            :π,                       # Built-in constant
            :Inf64,                   # Built-in constant
        ])
            @eval export $sym
        end
    end
end

