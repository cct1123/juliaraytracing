using LinearAlgebra

# Define the vector structure ===============================================================================
# mutable struct Vec3 <: AbstractVector{Float64}
#     data::NTuple{3, Float64}
# end
# # Overloading basic operators for Vec3
# Base.size(::Vec3) = (3,)
# Base.getindex(v::Vec3, i::Int) = v.data[i]
# Base.setindex!(v::Vec3, x, i::Int) = Vec3(Base.setindex(v.data, x, i))  # Immutable update
# # Add two vectors
# Base.:+(a::Vec3, b::Vec3) = Vec3(a.data .+ b.data)
# # Subtract two vectors
# Base.:-(a::Vec3, b::Vec3) = Vec3(a.data .- b.data)
# # Multiply vector by scalar
# Base.:*(s::Number, v::Vec3) = Vec3(s .* v.data)
# Base.:*(v::Vec3, s::Number) = Vec3(v.data .* s)
# # Divide vector by scalar
# /(v::Vec3, t::Number) = Vec3(v.data./t)
# # Define dot product
# dot(a::Vec3, b::Vec3) = sum(a.data .* b.data)
# ⋅(a::Vec3, b::Vec3) = dot(a, b)
# # Define cross product
# cross(a::Vec3, b::Vec3) = Vec3(
#     a.data[2] * b.data[3] - a.data[3] * b.data[2],
#     a.data[3] * b.data[1] - a.data[1] * b.data[3],
#     a.data[1] * b.data[2] - a.data[2] * b.data[1]
# )
# ×(a::Vec3, b::Vec3) = cross(a, b)
# # Define norm (magnitude)
# norm(v::Vec3) = sqrt(sum(v.data .^ 2))
# # Define magnitude (same as norm)
# magnitude(v::Vec3) = norm(v)
# # Define magnitude (same as norm)
# function normalize(v::Vec3)
#     mag = norm(v)
#     mag == 0.0 ? Vec3(0.0, 0.0, 0.0) : v / mag
# end
# Base.show(io::IO, v::Vec3) = print(io, "Vec3(", v.data[1], ", ", v.data[2], ", ", v.data[3], ")")

const Vec3 = Vector{Float64}
const Vec2 = Vector{Float64}
const Point = Vec3
# ========================================================================================================

# Define the Ray structure ===============================================================================
mutable struct Ray
    origin::Point # 3d vector origin position
    velocity::Vec3 # 3d vector velocity with magnitdue defined c=1
    polarization:: Vec3 # 2d vector polariation looking into the propagation direction

    function Ray(origin::Point, velocity::Vec3, polarization::Vec3)
        new(origin, velocity, )  # Normalize direction for consistency
    end

    # function Ray(origin::Point, velocity::Vec3, polarization::Vec3)
    #     new(origin, velocity, )  # Normalize direction for consistency
    # end
end
# ==============================================================================================================

# Define the structure for light sources ===============================================================================
struct Source
    # define the point source structure
    center::Point # 3d vector center position
    intensity::Float64 # intensity of the light source
    function Source(origin::Point, velocity::Vec3, polarization::Vec2)
        new(origin, velocity, )  # Normalize direction for consistency
    end
end

struct EnsembleSource <: Source
    objects::Vector{Source}
end

# ==============================================================================================================

# Define the structure for material, boundaries, and objects ===============================================================================
# Abstract material 
abstract type Material end

# Metal material with fuzziness
struct Metal <: Material
    # fake metal material
    albedo::Vec3  # Surface color (reflectance)
    fuzz::Float64 # Fuzziness factor (0 = no fuzz, 1 = max fuzz)
end

struct Dielectric <: Material
    # define simple dielectric material
    # TODO: to include complex dielectric instead of a simple refractive index
    n::Float64  # Refractive index of the material
end


# Abstract Object for the rays to hit 
abstract type AbtractObject end
abstract type AbstractBoundary end

struct  Frame 
    origin:: Point
    orientation:: Matrix{Float64} # a 3x3 rotation matrix that transform the lab frame to the object frame
    
end

# example shape function: Define a spherical surface as the shape function
function spherical_surface(center::Vec3, radius::Float64)
    return (p::Vec3) -> norm(p - center) - radius  # Sphere equation: ||p - center|| = radius
end

struct Boundary <: AbstractBoundary
    frame:: Frame
    n1:: Material # the material below the shape
    n2:: Material # the material above the shape
    shape:: Function
    bounds:: Tuple{}

end
struct BoundaryList <: Boundary
    boundaries: Vector{Boundary}
end

mutable struct HitRecord
end

struct Object <: AbtractObject
    frame:: Frame # the object frame


end


struct BoxObject <: Object
    center:: Point 
    size:: Vec3 
    material:: Material

end