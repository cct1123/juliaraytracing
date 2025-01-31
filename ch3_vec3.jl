mutable struct Vec3
    x::Float64
    y::Float64
    z::Float64

    function Vec3(x::Float64, y::Float64, z::Float64)
        new(x, y, z)
    end
end
const Point3 = Vec3 # alias Point3 to Vec3

# Overloading basic operators for Vec3
import Base: +, -, *, /, ⋅, ×

# Add two vectors
+(v1::Vec3, v2::Vec3) = Vec3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z)

# Subtract two vectors
-(v1::Vec3, v2::Vec3) = Vec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z)

# Multiply vector by scalar
*(v::Vec3, t::Float64) = Vec3(v.x * t, v.y * t, v.z * t)
*(t::Float64, v::Vec3) = Vec3(v.x * t, v.y * t, v.z * t)

# Divide vector by scalar
/(v::Vec3, t::Float64) = Vec3(v.x / t, v.y / t, v.z / t)

# Dot product
function dot(v1::Vec3, v2::Vec3)
    v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
end

⋅(v1::Vec3, v2::Vec3) = dot(v1, v2)

# Cross product
function cross(v1::Vec3, v2::Vec3)
    Vec3(
        v1.y * v2.z - v1.z * v2.y,
        v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x
    )
end
×(v1::Vec3, v2::Vec3) = cross(v1, v2)

# Magnitude of a vector
function norm(v::Vec3)
    sqrt(v ⋅ v)
end

# Normalize a vector
function normalize(v::Vec3)
    mag = norm(v)
    mag == 0.0 ? Vec3(0.0, 0.0, 0.0) : v / mag
end

# Print vector nicely
Base.show(io::IO, v::Vec3) = print(io, "Vec3(", v.x, ", ", v.y, ", ", v.z, ")")

# Example Usage
v1 = Vec3(3.0, 4.0, 0.0)
v2 = Vec3(1.0, 2.0, 3.0)
if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    println("v1 + v2 = ", v1 + v2)         # Vector addition
    println("v1 - v2 = ", v1 - v2)         # Vector subtraction
    println("Dot product = ", v1 ⋅ v2) # Dot product
    println("Cross product = ", v1 × v2)   # Cross product
    println("Magnitude of v1 = ", norm(v1)) # Magnitude
    println("Normalized v1 = ", normalize(v1)) # Normalized vector
end
