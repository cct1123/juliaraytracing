
module MathModule # module begin--------------------------

using LinearAlgebra
const Vec3 = Vector{Float64};
const Vec2 = Vector{Float64};
const Point = Vec3;
const Direction = Vec3;
const invalid_vector =[NaN, NaN, NaN]

struct  Frame 
    origin:: Point
    orientation:: Matrix{Float64} # a 3x3 rotation matrix that transform the lab frame to the object frame
    function Frame(origin::Point=[0.0, 0.0, 0.0], orientation::Matrix{Float64}=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
        new(origin, orientation)
    end
end

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


# how the rays interact with the objects and surfaces=============================================================================
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


function transform_direction_to(vec::Vector{Float64}, frame::Frame)
    return transpose(frame.orientation * vec)
end

function transform_direction_from(vec::Vector{Float64}, frame::Frame)
    return frame.orientation * vec
end

function transform_point_to(vec::Vector{Float64}, frame::Frame)
    return frame.orientation * (vec - frame.origin)
end

function transform_point_from(vec::Vector{Float64}, frame::Frame)
    return frame.orientation * vec + frame.origin
end

function reflect(normal::Vector, incident::Vector)
    cosI = -dot(normal, incident)
    return incident + 2 * cosI * normal
end

function refract(normal::Vector, incident::Vector, n1::Float64, n2::Float64)
    # normal is pointing from n2 to n1
    # incident ray is travelling from n1 to n2

    n = n1 / n2
    cosI = -dot(normal, incident)
    if cosI < 0
        n = 1/n
        cosI = -1.0*cosI
        normal = -1.0*normal 
    end
    sinT2 = n^2 * (1.0 - cosI^2)
    if sinT2 > 1.0  
        return invalid_vector  # TIR
    end
    cosT = sqrt(1.0 - sinT2)
    return n * incident + (n * cosI - cosT) * normal
end

function reflectance(normal::Vector, incident::Vector, n1::Float64, n2::Float64)
    # normal is pointing from n2 to n1
    # incident ray is travelling from n1 to n2
    n = n1 / n2
    cosI = -dot(normal, incident)
    if cosI < 0
        n = 1/n
        cosI *= -1.0
    end
    sinT2 = n^2 * (1.0 - cosI^2)
    if sinT2 > 1.0 
        return 1.0  # TIR
    end
    cosT = sqrt(1.0 - sinT2)
    rOrth = (n * cosI -  cosT) / (n * cosI + cosT)
    rPar = (cosI - n * cosT) / ( cosI + n * cosT)
    return (rOrth^2 + rPar^2) / 2.0
end

function fresnel(normal::Vector, incident::Vector, n1::Float64, n2::Float64)
    # normal is pointing from n2 to n1
    # incident ray is travelling from n1 to n2
    n = n1 / n2
    cosI = -dot(normal, incident)
    if cosI < 0
        n = 1/n
        cosI = -1.0*cosI
        normal = -1.0*normal 
    end
    sinT2 = n^2 * (1.0 - cosI^2)

    if sinT2 > 1.0 
        reflected = incident + 2 * cosI * normal
        # # WARNING the ray would be trapped forever if we don't consider absorption of material
        # return 0.0, normalize(incident+cosI*normal), 1.0, reflected  # TIR
        
        return 0.0, normalize(incident+cosI*normal), 1.0, reflected  # TIR
    end

    cosT = sqrt(1.0 - sinT2)
    # for unpolarized light
    rOrth = (n * cosI -  cosT) / (n * cosI + cosT)
    rPar = (cosI - n * cosT) / ( cosI + n * cosT)
    bigr = (rOrth^2 + rPar^2) / 2.0
    refracted = n * incident + (n * cosI - cosT) * normal
    reflected = incident + 2 * cosI * normal
    return 1.0-bigr, refracted, bigr, reflected
end

function rSchlick2(normal::Vector, incident::Vector, n1::Float64, n2::Float64)
    r0 = (n1 - n2) / (n1 + n2)
    r0 = r0 * r0  # Fix original code's typo: r0*r0 -> r0 = r0^2
    cosX = -dot(normal, incident)
    if n1 > n2
        n = n1 / n2
        sinT2 = n^2 * (1.0 - cosX^2)  # Fix original code's cosI -> cosX
        if sinT2 > 1.0
            return 1.0
        end
        cosX = sqrt(1.0 - sinT2)
    end
    x = 1.0 - cosX
    return r0 + (1.0 - r0) * x^6  # Original: x*x*x*x*x*x
end

# Define a function to compute ray-surface intersection
function find_intersection(f::Function, o::Point, d::Vector{Float64}, t_min=1e-6, t_max=1e3)
    """
    Finds the intersection of a ray with an implicit surface f(x, y, z) = 0
    using root-finding.

    Arguments:
    - f: function defining the implicit surface, f(x, y, z)
    - o: ray origin (3-element tuple or vector)
    - d: ray direction (3-element tuple or vector, should be normalized)
    - t_min, t_max: search range for t

    Returns:
    - t_intersect: the smallest positive intersection t (or nothing if no intersection)
    """
    # Search for a root in the range [t_min, t_max]
    try
        roots = find_zeros(t -> f(o + t * d), t_min, t_max)

        # Find the smallest root
        t_intersect = minimum(roots)

        # t_intersect = find_zero((t) -> f(o + t * d), (t_min, t_max), Roots.AlefeldPotraShi())  # Robust root-finding
        return t_intersect > 0 ? t_intersect : nothing  # Ensure positive intersection
    catch
        return nothing  # No valid root found
    end
end

export Point, Direction, Vec3, Vec2, Frame, rotation_matrix, rotate

end # module end-------------------------------


